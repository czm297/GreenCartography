#include "GreenCartographyManager.h"
#include "NSGA2/NSGA2.h"
#include "EnergyReduction/pictureOp.h"
#include "EnergyReduction/constants.h"
#include <sstream>
#include "SemanticExtract/MapRelationship.h"


GreenCartographyManager::GreenCartographyManager()
{
}

GreenCartographyManager::~GreenCartographyManager()
{
}


bool GreenCartographyManager::getMapfile()
{
	m_ctAllCategory2.clear();
	clusters.clear();
	IsCustoms.clear();
	fgValues.clear();

	std::ifstream inFile(mapfilepath, std::ios::in);
	std::string lineStr;
	ctCategory layer0;
	layer0.list.resize(1);
	ctColor c0;
	while (getline(inFile, lineStr)) 
	{
		std::stringstream ss(lineStr);
		std::string str;
		std::vector<std::string> lineArray;
		while (getline(ss, str, ','))
			lineArray.push_back(str);

		layer0.id = lineArray[0];
		layer0.type2 = atoi(lineArray[1].c_str());
		c0.r = atoi(lineArray[2].c_str());
		c0.g = atoi(lineArray[3].c_str());
		c0.b = atoi(lineArray[4].c_str());
		c0.proportion = atof(lineArray[5].c_str());
		c0.type2 = layer0.type2;

		layer0.list[0] = c0;
		layer0.proportion = c0.proportion;
		m_ctAllCategory2.emplace_back(layer0);
	}
	inFile.close();
	if (m_ctAllCategory2.size() == 0)
		return false;

	string relationshipFile;
	map<int, vector<int>> relationships; 
	PictureOP::getRelationshipPath(mapfilepath, relationshipFile); 
	std::ifstream inFile2(relationshipFile, std::ios::in);
	int index = 0;

	while (getline(inFile2, lineStr))
	{
		std::stringstream ss(lineStr);
		std::string str;
		std::vector<std::string> lineArray;
		while (getline(ss, str, ','))
			lineArray.push_back(str);
		if (atoi(lineArray[1].c_str()) == -2)
		{

		}
		else
		{
			int id = atoi(lineArray[1].c_str());
			map<int, vector<int>>::iterator iter = relationships.find(id);
			if (iter != relationships.end())
			{
				iter->second.emplace_back(index); 
			}
			else  
			{
				vector<int> list;
				list.emplace_back(index);
				relationships.insert(pair<int, vector<int>>(atoi(lineArray[1].c_str()), list));
			}
		}
		IsCustoms.emplace_back(atoi(lineArray[2].c_str()));
		fgValues.emplace_back(atof(lineArray[3].c_str()));

		index++;
	}
	inFile2.close();
	for (auto r : relationships)
		clusters.emplace_back(r.second);

	return true;
}

bool GreenCartographyManager::EnergyReduction2()
{
	if (m_edNSGAPopNum % 2 == 1) 
		m_edNSGAPopNum++;

	NSGA2* mnsga = new NSGA2();
	mnsga->mapName = PictureOP::getMapName(mapfilepath);
	mnsga->NSGA2popsize = m_edNSGAPopNum;
	mnsga->FeaLimitation = 5;
	mnsga->NSGA2generation = m_edIterationTime;
	mnsga->iniColorDis(mcolorSDis, mColorDis3, mColorDis4, mColorDis5,mColorDis6);
	mnsga->iniParameters(m_ctAllCategory2, clusters, IsCustoms, fgValues);

	population pop(mnsga);
	pop.maincal();
	resG.resize(m_ResultNum);
	resultScores.resize(m_ResultNum);
	for (int i = 0; i < m_ResultNum; i++)
	{
		resG[i] = pop.getResult(i * 1.0 / (m_ResultNum - 1), resultScores[i]); 
	}
	E = mnsga->iniEnergy;
	Output(pop);
	return true;
}

void GreenCartographyManager::Output(const population& pop)
{
	if (!PictureOP::exists("extractedFeature"))
	{
		int ret = _mkdir("extractedFeature"); 
		if (ret != 0)
			return;
	}
	for (int j = 0; j < m_ResultNum; j++)
	{
		string  outp = "extractedFeature/MapColor" + population::doubleToString((j * 1.0 / (m_ResultNum - 1))) + ".csv";
		auto& res = resG[j];
		cout << res.size();
		ofstream xxx(outp);
		for (int i = 0; i < m_ctAllCategory2.size(); i++)
		{
			xxx << m_ctAllCategory2[i].id << "," << m_ctAllCategory2[i].type2 << "," 
				<< m_ctAllCategory2[i].list[0].r << "," << m_ctAllCategory2[i].list[0].g << "," << m_ctAllCategory2[i].list[0].b 
				<< ","
				<< res[i].r << "," << res[i].g << "," << res[i].b << endl;
		}
		xxx.close();
	}

	FILE* p;
	p = fopen("extractedFeature/My_NSGA2.csv", "w+");
	for (int i = 0; i < m_edNSGAPopNum; i++)
	{		
		fprintf(p, "%f,%f,%f,%f,%f\n", 1 - pop.P[i].fvalue[0], pop.P[i].fvalue[1], pop.P[i].IsFeasible, pop.P[i].mNSGA->iniEnergy, 100 * (1 - pop.P[i].fvalue[1] / pop.P[i].mNSGA->iniEnergy));
	}
	fclose(p);
}

vector<vector<double>>GreenCartographyManager::getResultScores()
{
	return resultScores;
}

string GreenCartographyManager::getResultPathByIndex(int index)
{
	if (index < 0)
		index = 0;
	if (index >= m_ResultNum)
		index = m_ResultNum - 1;

	float rindex = index * 1.0 / (m_ResultNum - 1);
	return  "extractedFeature/result_" + population::doubleToString(rindex) + ".csv";
}


void GreenCartographyManager::MapSemanticExtract(const string& mapfile, float eps)
{
	MapRelationship mp;
	mp.setDBSCANEps(eps);
	std::vector<std::vector<double>> rgbs;
	std::vector<double> color(3);
	std::vector<string> layernames;
	std::vector<bool> isCus;
	std::vector<float> bgs;
	std::ifstream inFile(mapfile, std::ios::in);
	std::string lineStr;

	int type = -1;
	while (getline(inFile, lineStr)) 
	{
		std::stringstream ss(lineStr);
		std::string str;
		std::vector<std::string> lineArray;

		while (getline(ss, str, ','))
			lineArray.push_back(str);

		color[0] = atoi(lineArray[2].c_str());
		color[1] = atoi(lineArray[3].c_str());
		color[2] = atoi(lineArray[4].c_str());

		type = atoi(lineArray[1].c_str());
		if (type == 3)
		{
			bgs.emplace_back(fgpointL);
		}
		else if (type == 1)
		{
			bgs.emplace_back(fglineL);
		}
		else if (type == 0)
		{
			bgs.emplace_back(fgpolygonL);
		}
		else
			bgs.emplace_back(0.0f);

		layernames.emplace_back(lineArray[0]);
		isCus.emplace_back(NSGA2::IsConstLayer(lineArray[0]));
		rgbs.emplace_back(color);
	}
	inFile.close();
	mp.calCluster(rgbs);

	saveMapRelationshipfile(mapfile, layernames, mp.clustersAll, isCus, bgs);
}

void GreenCartographyManager::saveMapRelationshipfile(const string& mapfile, const std::vector<string>& layernames, const vector<vector<int>>& clustersAll,const std::vector<bool>& iscus, const std::vector<float>& bgs)
{
	string rpath;
	PictureOP::getRelationshipPath(mapfile, rpath);
	ofstream result(rpath);
	for (int i = 0; i < layernames.size(); i++)
	{
		result << layernames[i] << "," << clustersAll[i][0] << "," << iscus[i] << "," << bgs[i] << endl;
	}
	result.close();
}
