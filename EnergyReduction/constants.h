#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "../base.h"

#include <string>
#include <vector>

const int constWindowWidth       = 400;
const int constWindowHeight      = 400;

const int constNumCellsAcross = 10;
const int constNumCellsDown = 10;

const int constNumCellsAcrossLine = 7;
const int constNumCellsDownLine = 7;

const int constNumCellsAcrossPoint = 15;
const int constNumCellsDownPoint = 15;

const int constNumCellsAcrossBack = 3;

//number of weights each node must contain. One for each element of 
//the input vector. In this example it is 3 because a color is
//represented by its red, green and blue components. (RGB)
const int     constSizeOfInputVector   = 4;

//the number of epochs desired for the training
const int    constNumIterations       = 3000;

//the value of the learning rate at the start of training
const double constStartLearningRate = 0.3;

/*   uncomment the following if you'd like the SOM to classify randomly created training sets  */

//#define RANDOM_TRAINING_SETS

typedef struct MAPSTYLETRANSFORM_EXPORT ctColor
{
	std::string id = "";

	int r = -1;
	int g = -1;
	int b = -1;

	double proportion = 0.0;
	int count = 0;
	std::string type;
	int type2 = -1; 
	bool selected = false;
	int x = -1;
	int y = -1;
	float LineWidth = -1.0;
	float LineLength = 0.0;

	double H;
	double S;
	double V;

	double L;
	double A;
	double B;
};


typedef struct MAPSTYLETRANSFORM_EXPORT ctSelectInfo
{
	std::string id = "";
	int x = -1;
	int y = -1;
	ctColor color;
};


typedef struct MAPSTYLETRANSFORM_EXPORT ctCategory
{
	std::string id = "";
	std::vector<ctColor> list;
	
	double proportion = 0.0;
	int type2 = -1; 
};


typedef struct MAPSTYLETRANSFORM_EXPORT ctRange
{
public:
	ctRange()
		:minx_(1e38), miny_(1e38), minz_(1e38), maxx_(-1e38), maxy_(-1e38), maxz_(-1e38)
	{

	};

	void extend_box(const double& x, const double& y)
	{
		minx_ = x < minx_ ? x : minx_;
		miny_ = y < miny_ ? y : miny_;

		maxx_ = x > maxx_ ? x : maxx_;
		maxy_ = y > maxy_ ? y : maxy_;
	}

	std::vector<int> center()
	{
		return{ (int)((minx_ + maxx_ )/2.0),(int)((miny_ + maxy_) / 2.0) };
	}

private:
	double minx_;
	double maxx_;
	double miny_;
	double maxy_;
	double minz_;
	double maxz_;
};


#ifdef RANDOM_TRAINING_SETS

const int    constMaxNumTrainingSets  = 20;
const int    constMinNumTrainingSets  = 5;



#endif

#endif