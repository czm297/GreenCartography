#include "NSGA2.h"
#include"opencv2/core/core.hpp"


void NSGA2::iniColorDis(const float& mcolorSDis, const float& mColorDis3, const float& mColorDis4, const float& mColorDis5, const float& mColorDis6)
{
	AsscoDis = mcolorSDis;
	delta2_2 = mColorDis3; 
	delta22 = mColorDis4; 
	delta4 = mColorDis5; 
	delta1 = mColorDis6;
}

void NSGA2::iniParameters(const std::vector<ctCategory>& InformationEntropy, const vector<vector<int>>& clusters, const vector<bool>& IsCustoms, const vector<float>& fgValues)
{
	Dimension2 = InformationEntropy.size();
	Clusters = clusters;
	m_dPictureInformation = InformationEntropy;
	std::vector<vector<double>>().swap(sourceLabs);

	for (auto& MapColor : m_dPictureInformation)
	{
		sourceLabs.emplace_back(ColorSpaceTransfer::RGB2Lab(MapColor.list[0].r, MapColor.list[0].g, MapColor.list[0].b));
	}

	iniEnergy = 0.0;
	for (int i = 0; i < m_dPictureInformation.size(); i++)
	{
		iniEnergy += ((2.5 * m_dPictureInformation[i].list[0].r / 255) + (2.5 * m_dPictureInformation[i].list[0].g / 255) + (7.5 *m_dPictureInformation[i].list[0].b / 255)) * m_dPictureInformation[i].proportion;
	}

	vector<int>().swap(type2);
	backLayerIDs.clear();

	for (auto& layer : m_dPictureInformation)
	{
		type2.emplace_back(layer.list[0].type2);
		if (layer.list[0].type2 == 2)
		{
			backLayerIDs.emplace_back(type2.size() - 1);
		}
	}

	constLayer.clear(); 
	conlayernum = 0;
	isColorConst.resize(m_dPictureInformation.size(), false);
	std::vector<float>().swap(fgValues_c);
	for (int i = 0; i < m_dPictureInformation.size(); i++)
	{
		if (IsCustoms[i])
		{
			conlayernum++;
			isColorConst[i] = true;
			constLayer.emplace_back(i);
		}
		fgValues_c.emplace_back(fgValues[i]);
	}

	mark.resize(NSGA2popsize);//A

	ofstream relationships("extractedFeature/relationship.csv");
	colorRelationShips.resize(Dimension2, vector<int>(Dimension2, -1));
	for (int i = 0; i < Dimension2 - 1; i++)
	{
		for (int j = i + 1; j < Dimension2; j++)
		{
			colorRelationShips[i][j] = getColorRelationship(i, j);
			colorRelationShips[j][i] = colorRelationShips[i][j];
			relationships << i << "," << j << "," << m_dPictureInformation[i].id << "," << m_dPictureInformation[j].id << "," << colorRelationShips[i][j] << endl;
		}
	}
}

bool NSGA2::IsConstLayer(const string& name)
{
	if (name.find("Water") != std::string::npos)
		return true;
	if (name.find("Green") != std::string::npos)
		return true;
	if (name.find("water") != std::string::npos)
		return true;
	if (name.find("grass") != std::string::npos)
		return true;
	if (name.find("wood") != std::string::npos)
		return true;
	if (name.find("green") != std::string::npos)
		return true;
	if (name.find("forest") != std::string::npos)
		return true;
	if (name.find("river") != std::string::npos)
		return true;
	if (name.find("ocean") != std::string::npos)
		return true;
	if (name.find("lake") != std::string::npos)
		return true;
	if (name.find("0 - 200") != std::string::npos)
		return true;
	if (name.find("ºþ") != std::string::npos)
		return true;
	if (name.find("Ë®Ìå") != std::string::npos)
		return true;
	if (name.find("ÂÌµØ") != std::string::npos)
		return true;
	if (name.find("Forests") != std::string::npos)
		return true;
	if (name.find("Deserts") != std::string::npos)
		return true; 
	if (name.find("Temperate humid") != std::string::npos)
		return true;
	return false;
}

double NSGA2::rand_real(double low, double high)
{
	double h;
	h = (high - low) * URAND + low + 0.001;
	if (h >= high)
		h = high - 0.001;
	return h;
}

int NSGA2::rand_int(int low, int high)
{
	return int((high - low + 1) * URAND) + low;
}

int NSGA2::getGroupID(int index)
{
	int num = 0;
	for (auto& cluster : Clusters)
	{
		if (find(cluster.begin(), cluster.end(), index) != cluster.end()) 
		{
			return num;
		}
		num++; 
	}
	return -1;
}

int NSGA2::getColorRelationship(const int& index1, const int& index2)//1:association¡¢2:order¡¢-1:difference
{
	for (auto& cluster : Clusters)
	{
		if (find(cluster.begin(), cluster.end(), index1) != cluster.end()
			&& find(cluster.begin(), cluster.end(), index2) != cluster.end())
		{
			if (cluster.size() == 2)
				return 1;
			else
				return 2;
		}
	}
	return -1;
}


void individual::init(NSGA2 nsga)
{
	value.resize(mNSGA->Dimension2);
	std::vector<vector<double>> temp;
	for (auto& MapColor : nsga.m_dPictureInformation)
	{
		temp.emplace_back(ColorSpaceTransfer::RGB2HSV(MapColor.list[0].r, MapColor.list[0].g, MapColor.list[0].b));
	}

	for (int i = 0; i < mNSGA->Dimension2; i++)
	{
		value[i].H = temp[i][0] + mNSGA->rand_real(-5, 5);
		value[i].S = 100 - temp[i][1];
		value[i].V = 100 - temp[i][2];

		if (value[i].H < 0)
			value[i].H = 0;
		else if (value[i].H > 360)
			value[i].H = 360;

		if (value[i].S < 0)
			value[i].S = 0;
		else if (value[i].S > 100)
			value[i].S = 100;

		if (value[i].V < 0)
			value[i].V = 0;
		else if (value[i].V > 100)
			value[i].V = 100;
	}
}

void individual::f_count()
{
	vector<vector<double>> colorsLab(mNSGA->Dimension2);
	vector<vector<double>> colorsHSV(mNSGA->Dimension2);

	vector<double>temp;
	for (int i = 0; i < mNSGA->Dimension2; i++)
	{
		colorsHSV[i] = { value[i].H,value[i].S,value[i].V };
		colorsLab[i] = ColorSpaceTransfer::HSB2Lab(value[i].H, value[i].S, value[i].V);

	}

	fvalue[0] = getFunction1_2(colorsHSV, colorsLab); 
	fvalue[1] = getFunction2(colorsHSV); 
	IsFeasible = getFunction3(colorsLab, 100 * (1 - fvalue[1] / mNSGA->iniEnergy)); 
}

void individual::getOrderGroupsDH(const vector<vector<double>>& colorsHSV, vector<float>& colorsdH)
{
	for (auto& cluster : mNSGA->Clusters)
	{
		float dh = 0.0f;

		if (cluster.size() >= 3) 
		{
			for (int i = 0; i < cluster.size() - 1; i++)
			{
				for (int j = i + 1; j < cluster.size(); j++)
				{
					dh += ColorSpaceTransfer::getColorDH(colorsHSV[i][0], colorsHSV[j][0]) / (j - i);
				}
			}
			colorsdH.emplace_back(dh);
		}
		else
			colorsdH.emplace_back(0.0f);
	}
}

double individual::getFunction1_2(const vector<vector<double>>& colorsHSV, const vector<vector<double>>& colorsLab)
{
	vector<vector<float>> Semantics(mNSGA->Dimension2, vector<float>(mNSGA->Dimension2, 0.0));
	double Semantic = 0.0;    
	double sem = 0.0;
	double convention = 0.0; 
	double objective1 = 0.0;
	int backlayerNum = mNSGA->backLayerIDs.size();
	int i = 0;

	for ( i = 0; i < mNSGA->Dimension2 - 1; i++)
	{
		for (int j = i + 1; j < mNSGA->Dimension2; j++)
		{
			if (mNSGA->type2[i] == 2 || mNSGA->type2[j] == 2)
			{
				if (mNSGA->type2[i] == 2 && mNSGA->type2[j] == 2)
				{
					float dis = sqrt(pow(colorsLab[i][0] - colorsLab[j][0], 2) + pow(colorsLab[i][1] - colorsLab[j][1], 2) + pow(colorsLab[i][2] - colorsLab[j][2], 2));
					sem = min(1.0f, dis / mNSGA->delta1); 
					Semantics[i][j] = sem;
					Semantics[j][i] = sem;
				}
				else 
				{
					float sem1 = 0;
					float dis = sqrt(pow(colorsLab[i][0] - colorsLab[j][0], 2) + pow(colorsLab[i][1] - colorsLab[j][1], 2) + pow(colorsLab[i][2] - colorsLab[j][2], 2));
					if (mNSGA->type2[i] == 2)
						sem1 = min(1.0f, dis / mNSGA->fgValues_c[j]);
					else if (mNSGA->type2[j] == 2)
						sem1 = min(1.0f, dis / mNSGA->fgValues_c[i]);

					float sem2 = 0;
					float dL = abs(colorsLab[i][0] - colorsLab[j][0]);
					sem2 = min(1.0f, dL / mNSGA->delta4);

					sem = sem1 * sem2;
					Semantics[i][j] = sem;
					Semantics[j][i] = sem;

				}
			}
			else
			{
				sem = getColorSemanticObjective(i, j, colorsHSV, colorsLab);
				Semantics[i][j] = sem;
				Semantics[j][i] = sem;
			}
		}
		Semantic += accumulate(Semantics[i].begin(), Semantics[i].end(), 0.0f) / (mNSGA->Dimension2 - 1);
	}

	int temp_m = 0;
	for (int i = 0; i < mNSGA->Dimension2; i++)
	{
		double tempresult = getColorConObjective(i, colorsLab);
		if (tempresult == 0)
			temp_m++;
		else
			convention += tempresult;
	}

	objective1 = Semantic / mNSGA->Dimension2 * convention / (mNSGA->Dimension2 - temp_m); 
	return 1 - objective1; 
}

double individual::getFunction2(const vector<vector<double>>& colorsHSV)
{
	int RGB[3] = { 0 };
	double HSV[3] = { 0 };
	double energy = 0;

	const int k = mNSGA->Dimension2;
	int(*rgb_all)[3] = new int[k][3];
	for (int i = 0; i < mNSGA->Dimension2; i++)
	{
		HSV[0] = colorsHSV[i][0];
		HSV[1] = colorsHSV[i][1];
		HSV[2] = colorsHSV[i][2];
		ColorSpaceTransfer::HSB2RGB(HSV, RGB);

		rgb_all[i][0] = RGB[0];
		rgb_all[i][1] = RGB[1];
		rgb_all[i][2] = RGB[2]; 
		energy += (2.5 * rgb_all[i][0] / 255 + 2.5 * rgb_all[i][1] / 255 + 7.5 * rgb_all[i][2] / 255) * mNSGA->m_dPictureInformation[i].proportion; 

		if (mNSGA->m_dPictureInformation[i].type2 == 2)
		{
			back_energy = (2.5 * rgb_all[i][0] / 255 + 2.5 * rgb_all[i][1] / 255 + 7.5 * rgb_all[i][2] / 255) * mNSGA->m_dPictureInformation[i].proportion;
		}
	}
	front_energy = energy - back_energy;

	return energy;
}

double individual::getColorSemanticObjective(const int& i, const int& j, const vector<vector<double>>& colorsHSV, const vector<vector<double>>& colorsLab)
{
	double dist = ColorSpaceTransfer::calColorDist(colorsLab[i], colorsLab[j]);
	float dH = ColorSpaceTransfer::getColorDH(colorsHSV[i][0], colorsHSV[j][0]);
	int classNum = j - i;

	float H1, H2;
	float S1, S2;
	float V1, V2;

	if (mNSGA->colorRelationShips[i][j] == 1)
	{
		return min(1.0, mNSGA->AsscoDis / dist); 
	}
	else if (mNSGA->colorRelationShips[i][j] == 2)
	{
		return min(1.0f, mNSGA->LdH / abs(dH))
			* min(min(1.0, colorsHSV[i][1] / (classNum * mNSGA->LdV2 + colorsHSV[j][1])), ((classNum + 1) * mNSGA->LdV2 + colorsHSV[j][1]) / colorsHSV[i][1])
			* min(min(1.0, colorsHSV[j][2] / (classNum * mNSGA->LdV2 + colorsHSV[i][2])), ((classNum + 1) * mNSGA->LdV2 + colorsHSV[i][2]) / colorsHSV[j][2]);
	}
	else
	{
		return min(1.0, dist / mNSGA->AsscoDis); 
	}
}

double individual::getColorConObjective(const int& index, const vector<vector<double>>& colorsLab)
{
	if (!mNSGA->isColorConst[index])
		return 0.0;
	else
	{
		float dis = sqrt(pow(colorsLab[index][0] - mNSGA->sourceLabs[index][0], 2) + pow(colorsLab[index][1] - mNSGA->sourceLabs[index][1], 2) + pow(colorsLab[index][2] - mNSGA->sourceLabs[index][2], 2));
		return min(1.0f, mNSGA->delta22 / dis);
	}
}

double individual::getFunction3(const vector<vector<double>>& colorsLab,double f2=-1) 
{
	double semantic = 0.0;
	double colortransform = 0.0;
	double ext = 0.0;
	int backlayerNum = mNSGA->backLayerIDs.size();

	for (int i = 0; i < mNSGA->Dimension2 - 1; i++)
	{
		for (int j = i + 1; j < mNSGA->Dimension2; j++)
		{
			float dis = sqrt(pow(colorsLab[i][0] - colorsLab[j][0], 2) + pow(colorsLab[i][1] - colorsLab[j][1], 2) + pow(colorsLab[i][2] - colorsLab[j][2], 2));
			if (dis > mNSGA->delta1)
				semantic += 0.0;
			else if (dis < 5)
				semantic += INT_MAX;
			else 
				semantic += ((mNSGA->delta1 - dis));  
		}
	}

	for (int i = 0; i < mNSGA->Dimension2; i++)
	{
		ext += max(0.0, colorsLab[i][0] - mNSGA->sourceLabs[i][0]);
		vector<double> colorsRGB(3);
		colorsRGB = ColorSpaceTransfer::LabToRgb(colorsLab[i][0], colorsLab[i][1], colorsLab[i][2]);

		if (colorsRGB[0] >= mNSGA->m_dPictureInformation[i].list[0].r || colorsRGB[1] >= mNSGA->m_dPictureInformation[i].list[0].g || colorsRGB[2] >= mNSGA->m_dPictureInformation[i].list[0].b)
		{
			ext += max(0.0, colorsRGB[0] - mNSGA->m_dPictureInformation[i].list[0].r);
			ext += max(0.0, colorsRGB[1] - mNSGA->m_dPictureInformation[i].list[0].g);
			ext += max(0.0, colorsRGB[2] - mNSGA->m_dPictureInformation[i].list[0].b);
		}
	}
	return semantic + ext + colortransform;
}

void individual::calObjective(vector<double>& resultscore)
{
	vector<vector<double>> colorsLab(mNSGA->Dimension2);
	vector<vector<double>> colorsHSV(mNSGA->Dimension2);

	vector<double>temp;
	for (int i = 0; i < mNSGA->Dimension2; i++)
	{
		colorsHSV[i] = { value[i].H,value[i].S,value[i].V };
		colorsLab[i] = ColorSpaceTransfer::HSB2Lab(value[i].H, value[i].S, value[i].V);

	}
	vector<double> f1;
	vector<double> f2;

	fvalue[0] = getFunction1_2(colorsHSV, colorsLab, f1);
	fvalue[1] = getFunction2(colorsHSV, f2);
	resultscore[0] = 1 - fvalue[0];
	resultscore[1] = fvalue[1];    
	resultscore[2] = f1[0];        
	resultscore[3] = f1[1];        
	resultscore[4] = f1[2];        
	resultscore[5] = mNSGA->iniEnergy;                                    
	resultscore[6] = 100 * (1 - fvalue[1] / mNSGA->iniEnergy);            
}

double individual::getFunction1_2(const vector<vector<double>>& colorsHSV, const vector<vector<double>>& colorsLab, vector<double>& f1)
{
	f1.resize(3);
	vector<vector<float>> Semantics(mNSGA->Dimension2, vector<float>(mNSGA->Dimension2, 0.0));
	double Semantic = 0.0;
	double sem = 0.0;
	double convention = 0.0;
	double objective1 = 0.0;
	double hierarchy = 0.0;
	float p = 1.0f;
	int backlayerNum = mNSGA->backLayerIDs.size();

	int i = 0;
	for (i = 0; i < mNSGA->Dimension2 - 1; i++)
	{
		for (int j = i + 1; j < mNSGA->Dimension2; j++)
		{
			if (mNSGA->type2[i] == 2 || mNSGA->type2[j] == 2)
			{
				if (mNSGA->type2[i] == 2 && mNSGA->type2[j] == 2) 
				{
					float dis = sqrt(pow(colorsLab[i][0] - colorsLab[j][0], 2) + pow(colorsLab[i][1] - colorsLab[j][1], 2) + pow(colorsLab[i][2] - colorsLab[j][2], 2));
					sem = min(1.0f, dis / mNSGA->delta1);
					Semantics[i][j] = sem;
					Semantics[j][i] = sem;
				}
				else 
				{
					double sem1 = 0;
					float dis = sqrt(pow(colorsLab[i][0] - colorsLab[j][0], 2) + pow(colorsLab[i][1] - colorsLab[j][1], 2) + pow(colorsLab[i][2] - colorsLab[j][2], 2));
					if (mNSGA->type2[i] == 2)
						sem1 = min(1.0f, dis / mNSGA->fgValues_c[i]);
					else if (mNSGA->type2[j] == 2)
						sem1 = min(1.0f, dis / mNSGA->fgValues_c[j]);
					double sem2 = 0;
					float dL = abs(colorsLab[i][0] - colorsLab[j][0]);
					sem2 = min(1.0f, dL / mNSGA->delta4);

					sem = sem1 * sem2;
					Semantics[i][j] = sem;
					Semantics[j][i] = sem;

				}
			}
			else  
			{
				sem = getColorSemanticObjective(i, j, colorsHSV, colorsLab);
				Semantics[i][j] = sem;
				Semantics[j][i] = sem;
			}
		}
		f1[0] += accumulate(Semantics[i].begin(), Semantics[i].end(), 0.0f) / (mNSGA->Dimension2 - 1);

	}
	f1[0] += accumulate(Semantics[i].begin(), Semantics[i].end(), 0.0f) / (mNSGA->Dimension2 - 1);
	f1[0] = f1[0] / mNSGA->Dimension2;

	int temp_m = 0;
	for (int i = 0; i < mNSGA->Dimension2; i++)
	{
		double tempresult = getColorConObjective(i, colorsLab);
		if (tempresult == 0)
			temp_m++;
		else
			convention += tempresult;
	}
	f1[1] = convention / (mNSGA->Dimension2 - temp_m);
	return 1 - f1[0] * f1[1];
}

double individual::getFunction2(const vector<vector<double>>& colorsHSV, vector<double>& f2)
{
	f2.resize(1);
	double w = 0;
	int RGB[3] = { 0 };
	double HSV[3] = { 0 };
	double energy = 0;

	const int k = mNSGA->Dimension2;
	int(*rgb_all)[3] = new int[k][3];
	for (int i = 0; i < mNSGA->Dimension2; i++)
	{
		HSV[0] = colorsHSV[i][0];
		HSV[1] = colorsHSV[i][1];
		HSV[2] = colorsHSV[i][2];
		ColorSpaceTransfer::HSB2RGB(HSV, RGB);

		rgb_all[i][0] = RGB[0];
		rgb_all[i][1] = RGB[1];
		rgb_all[i][2] = RGB[2];

		energy += (2.5 * rgb_all[i][0] / 255 + 2.5 * rgb_all[i][1] / 255 + 7.5 * rgb_all[i][2] / 255) * mNSGA->m_dPictureInformation[i].proportion;
	}

	f2[0] = energy;

	return f2[0];
}


//NSGA2
population::population(NSGA2* nsga)
{
	len_f = 0;
	mNSGA = nsga;
	F = new individual * [2 * mNSGA->NSGA2popsize];
	for (int i = 0; i < 2 * mNSGA->NSGA2popsize; i++)
	{
		F[i] = new individual[2 * mNSGA->NSGA2popsize];	
	}

	P.resize(mNSGA->NSGA2popsize);
	Q.resize(mNSGA->NSGA2popsize);
	R.resize(2 * mNSGA->NSGA2popsize);
	len.resize(2 * mNSGA->NSGA2popsize);

	int i;
	for (i = 0; i < mNSGA->NSGA2popsize; i++)
	{
		P[i].mNSGA = mNSGA;
		P[i].init(*nsga); 
		P[i].sp.resize(2 * mNSGA->NSGA2popsize);
	}
	for (i = 0; i < mNSGA->NSGA2popsize; i++)
	{
		P[i].f_count();
	}
	for (int i = 0; i < mNSGA->NSGA2popsize; i++)
	{
		Q[i].value.resize(mNSGA->Dimension2);
		Q[i].sp.resize(2 * mNSGA->NSGA2popsize);
		R[i].value.resize(mNSGA->Dimension2);
		R[mNSGA->NSGA2popsize + i].value.resize(mNSGA->Dimension2);
		R[i].sp.resize(2 * mNSGA->NSGA2popsize);
		R[mNSGA->NSGA2popsize + i].sp.resize(2 * mNSGA->NSGA2popsize);
	}
	Pnum = mNSGA->NSGA2popsize;
	Qnum = 0;
	Rnum = 0;
}

void population::maincal()
{
	int s, i, j;

	s = mNSGA->NSGA2generation; 
	make_new_pop();
	ofstream it("extractedFeature/iteration.csv");
	while (s--)
	{
		set_p_q(); 

		fast_nondominated_sort();
		Pnum = 0;
		i = 0;

		while (Pnum + len[i] <= mNSGA->NSGA2popsize)
		{
			calu_crowding_distance(i);
			for (j = 0; j < len[i]; j++)
				P[Pnum++] = F[i][j];
			i++;
			if (i >= len_f)break;
		}

		if (i < len_f)
		{
			calu_crowding_distance(i);
			f_sort(i); 
		}
		for (j = 0; j < mNSGA->NSGA2popsize - Pnum; j++)
			P[Pnum++] = F[i][j];

		if (s % 10 == 0)
		{
			int f1index = getMinIndex(0);
			int f2index = getMinIndex(1);
				
			it << mNSGA->NSGA2generation - s << "," << 1-F[0][f1index].fvalue[0] << "," << F[0][f1index].IsFeasible 
				<< "," << F[0][f2index].fvalue[1] << "," << F[0][f2index].IsFeasible << endl;
		}
		make_new_pop(); 
	}
	it.close();
}

std::vector<ctColor> population::getResult(const float& Rindex, vector<double>& resultscore)
{
	int Pindex = 0;
	if (Rindex == 0) 
		qsort(&P[0], mNSGA->NSGA2popsize, sizeof(individual), cmp1);
	else if (Rindex == 1)
		qsort(&P[0], mNSGA->NSGA2popsize, sizeof(individual), cmp2);
	else
	{
		qsort(&F[0][0], len[0], sizeof(individual), cmp1); 
		for (int i = 0; i < mNSGA->NSGA2popsize; i++)
		{
			if (P[i].fvalue[0] == F[0][int(len[0] * Rindex)].fvalue[0] && P[i].fvalue[1] == F[0][int(len[0] * Rindex)].fvalue[1])
			{
				Pindex = i;
				break;
			}
		}
	}
	resultscore.resize(7);
	P[Pindex].calObjective(resultscore);

	std::vector<ctColor> res;
	int colorsRGB[3];
	double hsv_temp[3];
	
	ofstream iii("extractedFeature/result_" + doubleToString(Rindex) + ".csv");
	
	double per = 0.0;
	per = (mNSGA->iniEnergy - P[Pindex].fvalue[1]) / mNSGA->iniEnergy;
	per = per * 100;

	iii << 1 - P[Pindex].fvalue[0] << "," << P[Pindex].fvalue[1] << "," << P[Pindex].IsFeasible << "," << mNSGA->iniEnergy << "," << per << '%'
		<< resultscore[2] << "," << resultscore[3] << "," << resultscore[4] << endl;

	res.resize(mNSGA->Dimension2);
	for (int i = 0; i < mNSGA->Dimension2; i++)
	{
		hsv_temp[0] = P[Pindex].value[i].H;
		hsv_temp[1] = P[Pindex].value[i].S;
		hsv_temp[2] = P[Pindex].value[i].V;

		ColorSpaceTransfer::HSB2RGB(hsv_temp, colorsRGB);

		iii << mNSGA->m_dPictureInformation[i].id << "," << mNSGA->type2[i] << "," << colorsRGB[0] << "," << colorsRGB[1] << "," << colorsRGB[2] << "," << mNSGA->m_dPictureInformation[i].proportion << std::endl;

		res[i].r = colorsRGB[0];
		res[i].g = colorsRGB[1];
		res[i].b = colorsRGB[2];
	}
	iii.close();

	return res;
}

void population::make_new_pop()
{
	int i, j, x, y, t1, t2, t3, t4;
	double  u, b;

	mNSGA->mark.assign(mNSGA->NSGA2popsize, 0);
	mNSGA->temp1.assign(mNSGA->NSGA2popsize, 0);
	t3 = 0;
	while (t3 < mNSGA->NSGA2popsize / 2)
	{
		while (t1 = t2 = mNSGA->rand_int(0, mNSGA->NSGA2popsize - 1), mNSGA->mark[t1]);
		while (t1 == t2 || mNSGA->mark[t2])
		{
			t2 = mNSGA->rand_int(0, mNSGA->NSGA2popsize - 1);
		}
		t1 = choice(t1, t2);
		mNSGA->temp1[t3++] = t1;
		mNSGA->mark[t1] = 1;
	}

	for (i = 0; i < mNSGA->NSGA2popsize; i += 2) 
	{
		for (j = 0; j < mNSGA->Dimension2; j++)
		{
			u = mNSGA->rand_real((0.0 + 1e-6), (1.0 - 1e-6));
			if (u <= 0.5)
				b = pow(2 * u, 1.0 / 2);
			else
				b = 1.0 / pow(2 * (1 - u), 1.0 / 2);

			x = y = mNSGA->rand_int(0, mNSGA->NSGA2popsize / 2 - 1);
			while (x == y)
				y = mNSGA->rand_int(0, mNSGA->NSGA2popsize / 2 - 1);

			Q[i].value[j].H = 1.0 / 2 * ((1 - b) * P[mNSGA->temp1[x]].value[j].H + (1 + b) * P[mNSGA->temp1[y]].value[j].H);
			Q[i].value[j].S = 1.0 / 2 * ((1 - b) * P[mNSGA->temp1[x]].value[j].S + (1 + b) * P[mNSGA->temp1[y]].value[j].S);
			Q[i].value[j].V = 1.0 / 2 * ((1 - b) * P[mNSGA->temp1[x]].value[j].V + (1 + b) * P[mNSGA->temp1[y]].value[j].V);

			if (Q[i].value[j].H < 0) 
				Q[i].value[j].H = 0;
			else if (Q[i].value[j].H > 360)
				Q[i].value[j].H = 360;

			if (Q[i].value[j].S < 0)
				Q[i].value[j].S = 0;
			else if (Q[i].value[j].S > 100)
				Q[i].value[j].S = 100;

			if (Q[i].value[j].V < 0)
				Q[i].value[j].V = 0;
			else if (Q[i].value[j].V > 100)
				Q[i].value[j].V = 100;


			u = mNSGA->rand_real((0.0 + 1e-6), (1.0 - 1e-6));
			if (u <= 0.5)
				b = pow(2 * u, 1.0 / 2);
			else
				b = 1.0 / pow(2 * (1 - u), 1.0 / 2);
			if (i + 1 < mNSGA->NSGA2popsize)
			{
				Q[i + 1].value[j].H = 1.0 / 2 * ((1 + b) * P[mNSGA->temp1[x]].value[j].H + (1 - b) * P[mNSGA->temp1[y]].value[j].H);
				Q[i + 1].value[j].S = 1.0 / 2 * ((1 + b) * P[mNSGA->temp1[x]].value[j].S + (1 - b) * P[mNSGA->temp1[y]].value[j].S);
				Q[i + 1].value[j].V = 1.0 / 2 * ((1 + b) * P[mNSGA->temp1[x]].value[j].V + (1 - b) * P[mNSGA->temp1[y]].value[j].V);

				if (Q[i + 1].value[j].H < 0)
					Q[i + 1].value[j].H = 0;
				else if (Q[i + 1].value[j].H > 360)
					Q[i + 1].value[j].H = 360;

				if (Q[i + 1].value[j].S < 0)
					Q[i + 1].value[j].S = 0;
				else if (Q[i + 1].value[j].S > 100)
					Q[i + 1].value[j].S = 100;

				if (Q[i + 1].value[j].V < 0)
					Q[i + 1].value[j].V = 0;
				else if (Q[i + 1].value[j].V > 100)
					Q[i + 1].value[j].V = 100;
			}
		}

		for (j = 0; j < mNSGA->Dimension2; j++)
		{
			int umH = 25;
			int umS = 10;
			int umV = 10;

			u = mNSGA->rand_real(0.0 + (1e-6), 1.0 - (1e-6));
			if (u < 0.5)
				u = pow(2 * u, 1.0 / (1 + umH)) - 1;
			else
				u = 1 - pow(2 * (1 - u), 1.0 / (1 + umH));
			Q[i].value[j].H = Q[i].value[j].H + u*360;


			u = mNSGA->rand_real(0.0 + (1e-6), 1.0 - (1e-6));
			if (u < 0.5)
				u = pow(2 * u, 1.0 / (1 + umS)) - 1;
			else
				u = 1 - pow(2 * (1 - u), 1.0 / (1 + umS));
			//if(u<0)
				Q[i].value[j].S = Q[i].value[j].S - u * 100;

			u = mNSGA->rand_real(0.0 + (1e-6), 1.0 - (1e-6));
			if (u < 0.5)
				u = pow(2 * u, 1.0 / (1 + umV)) - 1;
			else
				u = 1 - pow(2 * (1 - u), 1.0 / (1 + umV));
			//if(u<0)
				Q[i].value[j].V = Q[i].value[j].V + u*100;

			if (Q[i].value[j].H < 0)
				Q[i].value[j].H = 0;
			else if (Q[i].value[j].H > 360)
				Q[i].value[j].H = 360;
			if (Q[i].value[j].S < 0)
				Q[i].value[j].S = 0;
			else if (Q[i].value[j].S > 100)
				Q[i].value[j].S = 100;
			if (Q[i].value[j].V < 0)
				Q[i].value[j].V = 0;
			else if (Q[i].value[j].V > 100)
				Q[i].value[j].V = 100;


			u = mNSGA->rand_real(0.0 + (1e-6), 1.0 - (1e-6));
			if (u < 0.5)
				u = pow(2 * u, 1.0 / (1 + umH)) - 1;
			else
				u = 1 - pow(2 * (1 - u), 1.0 / (1 + umH));
			Q[i + 1].value[j].H = Q[i + 1].value[j].H + u*360;


			u = mNSGA->rand_real(0.0 + (1e-6), 1.0 - (1e-6));
			if (u < 0.5)
				u = pow(2 * u, 1.0 / (1 + umS)) - 1;
			else
				u = 1 - pow(2 * (1 - u), 1.0 / (1 + umS));
			//if(u<0)
				Q[i + 1].value[j].S = Q[i + 1].value[j].S - u * 100;

			u = mNSGA->rand_real(0.0 + (1e-6), 1.0 - (1e-6));
			if (u < 0.5)
				u = pow(2 * u, 1.0 / (1 + umV)) - 1;
			else
				u = 1 - pow(2 * (1 - u), 1.0 / (1 + umV));
			//if(u<0)
				Q[i + 1].value[j].V = Q[i + 1].value[j].V + u*100;

			if (Q[i + 1].value[j].H < 0)
				Q[i + 1].value[j].H = 0;
			else if (Q[i + 1].value[j].H > 360)
				Q[i + 1].value[j].H = 360;
			if (Q[i + 1].value[j].S < 0)
				Q[i + 1].value[j].S = 0;
			else if (Q[i + 1].value[j].S > 100)
				Q[i + 1].value[j].S = 100;
			if (Q[i + 1].value[j].V < 0)
				Q[i + 1].value[j].V = 0;
			else if (Q[i + 1].value[j].V > 100)
				Q[i + 1].value[j].V = 100;
		}
	}

	Qnum = mNSGA->NSGA2popsize;
	for (i = 0; i < mNSGA->NSGA2popsize; i++)
	{
		Q[i].mNSGA = mNSGA;
	}
	for (i = 0; i < mNSGA->NSGA2popsize; i++)
	{
		Q[i].f_count();
	}

}

int population::choice(int a, int b)
{
	if (P[a].rank < P[b].rank)
		return a;
	else if (P[a].rank == P[b].rank)
	{
		if (P[a].crowding_distance > P[b].crowding_distance)
			return a;
		else
			return b;
	}
	else
		return b;
}

void population::set_p_q() 
{
	Rnum = 0;
	Qnum = mNSGA->NSGA2popsize;
	int i;
	for (i = 0; i < Pnum; i++)
		R[Rnum++] = P[i];
	for (i = 0; i < Qnum; i++)
		R[Rnum++] = Q[i];
}

void population::fast_nondominated_sort()
{
	int i, j, k;
	vector<individual> H(2 * mNSGA->NSGA2popsize);

	int h_len = 0;
	for (i = 0; i < 2 * mNSGA->NSGA2popsize; i++)
	{
		R[i].np = 0; 
		R[i].is_dominated = 0; 
		len[i] = 0;
	}
	for (i = 0; i < 2 * mNSGA->NSGA2popsize; i++)
	{
		for (j = 0; j < 2 * mNSGA->NSGA2popsize; j++)
		{
			if (i != j)
			{
				if (e_is_dominated(R[i], R[j]))
					R[i].sp[R[i].is_dominated++] = j;
				else if (e_is_dominated(R[j], R[i]))
					R[i].np += 1;
			}
		}
		if (R[i].np == 0)
		{
			len_f = 1;
			F[0][len[0]++] = R[i];
		}

	}
	i = 0;
	while (len[i] != 0)
	{
		h_len = 0;
		for (j = 0; j < len[i]; j++)
		{
			for (k = 0; k < F[i][j].is_dominated; k++)
			{
				R[F[i][j].sp[k]].np--;
				if (R[F[i][j].sp[k]].np == 0)
				{
					H[h_len++] = R[F[i][j].sp[k]];
					R[F[i][j].sp[k]].rank = i + 2;
				}
			}
		}
		++i;
		if (i >= 2 * mNSGA->NSGA2popsize)
			break;
		len[i] = h_len;
		if (h_len != 0)
		{
			len_f++;
			for (j = 0; j < len[i]; j++)
				F[i][j] = H[j];
		}
	}
}

bool population::e_is_dominated(const individual& a, const individual& b)
{
	if (a.IsFeasible > mNSGA->FeaLimitation || b.IsFeasible > mNSGA->FeaLimitation)
	{
		if (a.IsFeasible < b.IsFeasible)
			return true;
		if (b.IsFeasible < a.IsFeasible)
			return false;
	}

	if ((a.fvalue[0] <= b.fvalue[0]) && (a.fvalue[1] <= b.fvalue[1])) 
	{
		if ((a.fvalue[0] == b.fvalue[0]) && a.fvalue[1] == b.fvalue[1])
			return false;
		else
			return true;
	}
	else
		return false;
}

void population::calu_crowding_distance(int i)
{
	int n = len[i];
	double m_max, m_min;
	int j;
	for (j = 0; j < n; j++)
		F[i][j].crowding_distance = 0;

	qsort(F[i], n, sizeof(individual), cmp1);
	F[i][0].crowding_distance = F[i][n - 1].crowding_distance = 0xffffff;
	m_max = -0xfffff;
	m_min = 0xfffff;
	for (j = 0; j < n; j++)
	{
		if (m_max < F[i][j].fvalue[0])
			m_max = F[i][j].fvalue[0];
		if (m_min > F[i][j].fvalue[0])
			m_min = F[i][j].fvalue[0];
	}
	for (j = 1; j < n - 1; j++)
		F[i][j].crowding_distance += (F[i][j + 1].fvalue[0] - F[i][j - 1].fvalue[0]) / (m_max - m_min);

	qsort(F[i], n, sizeof(individual), cmp2);
	F[i][0].crowding_distance = F[i][n - 1].crowding_distance = 0xffffff;
	m_max = -0xfffff;
	m_min = 0xfffff;
	for (j = 0; j < n; j++)
	{
		if (m_max < F[i][j].fvalue[1])
			m_max = F[i][j].fvalue[1];
		if (m_min > F[i][j].fvalue[1])
			m_min = F[i][j].fvalue[1];
	}
	for (j = 1; j < n - 1; j++)
		F[i][j].crowding_distance += (F[i][j + 1].fvalue[1] - F[i][j - 1].fvalue[1]) / (m_max - m_min);
}

int population::cmp1(const void* a, const void* b)
{
	const individual* e = (const individual*)a;
	const individual* f = (const individual*)b;
	if (e->fvalue[0] == f->fvalue[0])
		return 0;
	else if (e->fvalue[0] < f->fvalue[0])
		return -1;
	else return 1;
}

int population::cmp2(const void* a, const void* b)
{
	const individual* e = (const individual*)a;
	const individual* f = (const individual*)b;
	if (e->fvalue[1] == f->fvalue[1])
		return 0;
	else if (e->fvalue[1] < f->fvalue[1])
		return -1;
	else return 1;
}

void population::f_sort(int i)
{
	int n;
	n = len[i];
	qsort(F[i], n, sizeof(individual), cmp_c_d);
}

int population::cmp_c_d(const void* a, const void* b)
{
	const individual* e = (const individual*)a;
	const individual* f = (const individual*)b;
	if (e->crowding_distance == f->crowding_distance)
		return 0;
	else if (e->crowding_distance < f->crowding_distance)
		return 1;
	else
		return -1;
}

int population::getMinIndex(const int& findex)
{
	double minF = F[0][0].fvalue[findex];
	int index = 0;
	for (int i = 1; i < len[0]; i++)
	{
		if (F[0][i].fvalue[findex] < minF && F[0][i].IsFeasible< mNSGA->FeaLimitation)
		{
			minF = F[0][i].fvalue[findex];
			index = i;
		}
	}
	return index;
}