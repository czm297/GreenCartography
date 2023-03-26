#pragma once

#include "NSGA2\NSGA2.h"
#include "EnergyReduction/constants.h"
#include <string>
#include <vector>

using namespace std;


class MAPSTYLETRANSFORM_EXPORT GreenCartographyManager
{
public:
	string mapfilepath;
	int m_edIterationTime = 2000; 
	int m_edNSGAPopNum = 200; 
	int m_ResultNum = 5;
	float mcolorSDis = 40;
	float mColorDis3 = 53;
	float mColorDis4 = 20;
	float mColorDis5 = 20;
	float mColorDis6 = 20;

	double E = 0.0;
	std::vector<ctCategory> m_ctAllCategory2;  
	std::vector<std::vector<int>> clusters;
	std::vector<bool> IsCustoms; 
	std::vector<float> fgValues; 

	float fgpointL = 32.0;
	float fglineL = 17.0f;
	float fgpolygonL = 10.5f;

	std::vector<std::vector<ctColor>> resG; 
	vector<vector<double>> resultScores;

	bool getMapfile();/
	bool EnergyReduction2();
	void Output(const population& pop);

	vector<vector<double>>  getResultScores(); 
	string getResultPathByIndex(int index);

	void MapSemanticExtract(const string& mapfile, float eps = 15.0f);
	void saveMapRelationshipfile(const string& mapfile, const std::vector<string>& layernames, const vector<vector<int>>& clustersAll, const std::vector<bool>& iscus, const std::vector<float>& bgs);

	GreenCartographyManager();
	~GreenCartographyManager();
};