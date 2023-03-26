#pragma once

#include <algorithm>
#include <vector>

#include "../EnergyReduction/constants.h"
#include "../EnergyReduction/ColorSapceTransfer.h"

#include <vector> 
#include <iostream>
#include <iostream>                                                                                                                                
#include <fstream>
#include <math.h>
#include <algorithm>
#include <omp.h>
#include <fstream>
#include "../base.h"


using namespace std;

#define URAND (rand()/(RAND_MAX+1.0))

struct colorHSV
{
	double H;
	double S;
	double V;
};

class NSGA2
{
public:
	string mapName;
	bool isBrightLimit = true;
	int NSGA2popsize; 
	int NSGA2generation;
	int FeaLimitation;
	double iniEnergy; 

	vector<vector<int>> Clusters;
	vector<ctCategory> m_dPictureInformation;
	int Dimension2;
	int conlayernum;
	vector<vector<double>> sourceLabs; 
	vector<int> mark;
	vector<int> temp1;
	vector<int> constLayer;
	vector<bool> isColorConst;  
	vector<vector<int>> colorRelationShips; 

	vector<int> type2; 
	vector<int> backLayerIDs;
	vector<float> fgValues_c;

	float AsscoDis = 40;
	float delta1 = 20;  
	float delta2 = 20;  
	float delta2_2 = 53; 
	float delta22 = 20;  
	float delta4 = 20;   

	float LdH = 9;     
	float LdH2 = 10;   
	float LdS = 10;
	float LdV = 8;
	float LdV2 = 10;

	void iniColorDis(const float& mcolorSDis, const float& mColorDis3, const float& mColorDis4, const float& mColorDis5, const float& mColorDis6); 
	void iniParameters(const std::vector<ctCategory>& InformationEntropy, const vector<vector<int>>& clusters, const vector<bool>& IsCustoms, const vector<float>& fgValues);

	double rand_real(double low, double high);
	int rand_int(int low, int high);
	int getColorRelationship(const int& index1, const int& index2); 
	int getGroupID(int index);
	static bool IsConstLayer(const string& name); 
};

class individual
{
private:
	vector<float> GroupdHs;
	double getColorSemanticObjective(const int& i, const int& j, const vector<vector<double>>& colorsHSV, const vector<vector<double>>& colorsLab);
	double getColorConObjective(const int& index, const vector<vector<double>>& colorsLab); 
	double getFunction1_2(const vector<vector<double>>& colorsHSV, const vector<vector<double>>& colorsLab)
	double getFunction1_2(const vector<vector<double>>& colorsHSV, const vector<vector<double>>& colorsLab, vector<double>& f1);
	double getFunction2(const vector<vector<double>>& colorsHSV);
	double getFunction2(const vector<vector<double>>& colorsHSV, vector<double>& f2);
	double getFunction3(const vector<vector<double>>& colorsLab, double f2); 
	void getOrderGroupsDH(const vector<vector<double>>& colorsHSV, vector<float>& colorsdH);

public:
	NSGA2* mNSGA;
	vector<int> sp;
	std::vector<colorHSV> value;
	int np;
	int is_dominated;
	double fvalue[2];       
	double IsFeasible = 0.0;
	double back_energy = 0.0;
	double front_energy = 0.0;
	int rank;
	double crowding_distance;
	void init(NSGA2);
	void f_count();
	void calObjective(vector<double>& resultscore);
};

class population
{
public:

	NSGA2* mNSGA;
	individual** F;
	vector<individual> P;
	vector<individual> Q;
	vector<individual> R;
	vector<int> len;
	int Rnum;
	int Pnum;
	int Qnum;
	int len_f;
	population(NSGA2*);
	void maincal();
	void make_new_pop();
	int choice(int a, int b);
	void set_p_q();
	void fast_nondominated_sort();
	void calu_crowding_distance(int i);
	void f_sort(int i);
	bool e_is_dominated(const individual& a, const individual& b);
	std::vector<ctColor> getResult(const float& index, vector<double>& resultscore);
	static string doubleToString(const double& dbNum)
	{
		char* chCode;
		chCode = new(std::nothrow)char[20];
		sprintf(chCode, "%.2lf", dbNum); 
		string strCode(chCode);
		delete[]chCode;
		return strCode;
	}

private:
	static int cmp1(const void *a, const void *b); 
	static int cmp2(const void *a, const void *b); 
	static int cmp_c_d(const void* a, const void* b);
	int getMinIndex(const int& findex);
};