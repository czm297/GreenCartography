#pragma once
#include <vector>
#include <string>
#include "../EnergyReduction/ColorSapceTransfer.h"
#include "../base.h"


struct metacolor
{
	int colorID;
	int clusterID;
	int clusterID2;
	bool selected = false;
	string name;
};

class DBSCANpoint {
public:


	int id;
	float z;
	float x;
	float y;
	int cluster = 0;

	int pointType = 1;//1 noise 2 border 3 core¡¢   0-diff 1-order 2-Asso 

	int pts = 0;//points in MinPts 
	int j = 0;
	std::vector<int> corepts;
	int visited = 0;
	DBSCANpoint() {}
	DBSCANpoint(float a, float b, float c, int d)
	{
		z = a;
		x = b;
		y = c;
		cluster = d;

	}
	double getDis(const DBSCANpoint & ot) {
		return sqrt((x - ot.x)*(x - ot.x) + (y - ot.y)*(y - ot.y) + (z - ot.z)*(z - ot.z));
	}
};
const int NOISE = -2;
const int NOT_CLASSIFIED = -1;

class DBSCANF
{
public:
	vector<DBSCANpoint> points;
private:
	int minPts;
	double eps;
	int size;
	vector<vector<int>> adjPoints;
	vector<bool> visited;
	vector<vector<int> > cluster;
	int clusterIdx;

	void dfs(int now, int c) {
		points[now].cluster = c;
		if (!isCoreObject(now)) return;

		for (auto&next : adjPoints[now]) {
			if (points[next].cluster != NOT_CLASSIFIED) continue;
			dfs(next, c);
		}
	}

	void checkNearPoints() {
		for (int i = 0; i<size; i++) {
			for (int j = 0; j<size; j++) {
				if (i == j) continue;
				if (points[i].getDis(points[j]) <= eps) {
					points[i].pts++;
					adjPoints[i].push_back(j);
				}
			}
		}
	}

	bool isCoreObject(int idx) {
		return points[idx].pts >= minPts;
	}

public:
	DBSCANF(double eps, int minPts, int clusterIdx, vector<DBSCANpoint> points) 
	{
		this->eps = eps;
		this->minPts = minPts;
		this->points = points;
		this->size = (int)points.size();
		adjPoints.resize(size);
		this->clusterIdx = clusterIdx;
	}

	void run() {
		checkNearPoints();
		for (int i = 0; i < size; i++) {
			points[i].cluster = -1;
		}
		for (int i = 0; i<size; i++) {
			if (points[i].cluster != NOT_CLASSIFIED) continue;

			if (isCoreObject(i)) {
				dfs(i, ++clusterIdx);
			}
			else {
				points[i].cluster = NOISE;
			}
		}

		cluster.resize(clusterIdx + 1);
		for (int i = 0; i<size; i++) {
			if (points[i].cluster != NOISE) {
				cluster[points[i].cluster].push_back(i);
			}
		}
	}

	vector<vector<int>> getCluster() {
		return cluster;
	}

	int getclusterIdx()
	{
		return clusterIdx;
	}
};

class MapRelationship
{
public:
	MapRelationship();
	~MapRelationship();

	vector<vector<int>> clustersAll;

	void calCluster(std::vector<std::vector<double>>& rgbs)
	{
		int D = rgbs.size();
		std::vector<DBSCANpoint> m_lab;
		inputColor.resize(D);
		for (int i = 0; i < D; i++)
		{
			DBSCANpoint temp;

			auto lab = ColorSpaceTransfer::RGB2Lab(rgbs[i][0], rgbs[i][1], rgbs[i][2]);


			inputColor[i].colorID= temp.id = i;
			temp.z = lab[0];
			temp.x = lab[1];
			temp.y = lab[2];
			
			m_lab.push_back(temp);
		}

		DBSCAN3(m_lab, DBSCANEps, DBSCANMinPts, DBSCANclustertime);

		map<int, int> clusters;

		for (int i = 0; i < D; i++)
		{
			cluster cc;
			

			if (m_lab[i].cluster == 0)
			{
				cc.type = 0;
				cc.colorIDs.emplace_back(i);
				colorCluster.emplace_back(cc);
			}
			else
			{
				map<int, int>::iterator it = clusters.find(m_lab[i].cluster);
				if (it != clusters.end())
				{
					colorCluster[(*it).second].colorIDs.emplace_back(i);
				}
				else
				{
					clusters[m_lab[i].cluster] = colorCluster.size();
					cc.colorIDs.emplace_back(i);
					colorCluster.emplace_back(cc);
				}
			}
		}

		for (int i = 0; i < colorCluster.size(); i++)
		{
			if (colorCluster[i].colorIDs.size() > 2)
				colorCluster[i].type = 1;
			if (colorCluster[i].colorIDs.size() == 2)
				colorCluster[i].type = 2;
		}
	}

	void setDBSCANEps(float eps)
	{
		DBSCANEps = eps;
	}

	float getDBSCANEps()
	{
		return DBSCANEps;
	}

private:

	float DBSCANEps = 25.0f;
	float DBSCANMinPts = 1;
	float DBSCANclustertime = 1;

	struct cluster
	{
		int type;  	
		vector<int> colorIDs;
	};

	vector<cluster> colorCluster;
	std::vector<metacolor> inputColor;

	void DBSCAN3(std::vector<DBSCANpoint>& datas, float Eps, int MinPts, int clustertime)
	{
		int clusterIdx = -1;
		clustersAll.resize(datas.size());


		DBSCANF d1(Eps, MinPts, clusterIdx, datas);
		d1.run();
		vector<vector<int>>  cluster0 = d1.getCluster();
		clusterIdx = d1.getclusterIdx();

		for (int i = 0; i < datas.size(); i++)
			inputColor[i].clusterID2 = d1.points[i].cluster;

		for (int i = 0; i < datas.size(); i++)
			clustersAll[i].emplace_back(d1.points[i].cluster);

		if (clustertime > 1)
		{
			for (int i = 0; i < cluster0.size(); i++)
			{
				if (cluster0[i].size()>2)
				{
					std::vector<DBSCANpoint> data2;
					for (int j = 0; j < cluster0[i].size(); j++)
						data2.emplace_back(datas[cluster0[i][j]]);

					DBSCANF d2(15, MinPts, -1, data2);
					d2.run();
					vector<vector<int>>  cluster1 = d2.getCluster();

					for (int j = 0; j < cluster0[i].size(); j++)
					{
						int index = cluster0[i][j];
						if (d2.points[j].cluster != -2)
						{
							inputColor[index].clusterID2 = d2.points[j].cluster;
						}
						clustersAll[index].emplace_back(d2.points[j].cluster);
					}

				}
			}
		}
	}

};



MapRelationship::MapRelationship()
{
}

MapRelationship::~MapRelationship()
{
}