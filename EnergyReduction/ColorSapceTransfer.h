#pragma once

#include <vector>
#include <list>
#include <algorithm>
#include <numeric>  
#include <tuple>
#include <map>
#include <math.h>

#ifndef _max
#define _max(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef _min
#define _min(a,b)            (((a) < (b)) ? (a) : (b))
#endif


using namespace std;

class MAPSTYLETRANSFORM_EXPORT ColorSpaceTransfer
{
public:
	static std::tuple<double, double, double> RGB2Lab1(double R, double G, double B)
	{
		double Lab_L, Lab_a, Lab_b;
		Lab_L = Lab_a = Lab_b = 0;

		double X, Y, Z;
		double r = R / 255.000; // rgb range: 0 ~ 1
		double g = G / 255.000;
		double b = B / 255.000;
		// gamma 2.2
		if (r > 0.04045)
			r = pow((r + 0.055) / 1.055, 2.4);
		else
			r = r / 12.92;
		if (g > 0.04045)
			g = pow((g + 0.055) / 1.055, 2.4);
		else
			g = g / 12.92;
		if (b > 0.04045)
			b = pow((b + 0.055) / 1.055, 2.4);
		else
			b = b / 12.92;
		// sRGB
		X = r * 0.436052025 + g * 0.385081593 + b * 0.143087414;
		Y = r * 0.222491598 + g * 0.716886060 + b * 0.060621486;
		Z = r * 0.013929122 + g * 0.097097002 + b * 0.714185470;
		// XYZ range: 0~100
		X = X * 100.000;
		Y = Y * 100.000;
		Z = Z * 100.000;
		// Reference White Point
		double ref_X = 96.4221;
		double ref_Y = 100.000;
		double ref_Z = 82.5211;
		X = X / ref_X;
		Y = Y / ref_Y;
		Z = Z / ref_Z;
		// Lab
		if (X > 0.008856)
			X = pow(X, 1 / 3.000);
		else
			X = (7.787 * X) + (16 / 116.000);
		if (Y > 0.008856)
			Y = pow(Y, 1 / 3.000);
		else
			Y = (7.787 * Y) + (16 / 116.000);
		if (Z > 0.008856)
			Z = pow(Z, 1 / 3.000);
		else
			Z = (7.787 * Z) + (16 / 116.000);
		Lab_L = (116.000 * Y) - 16.000;
		Lab_a = 500.000 * (X - Y);
		Lab_b = 200.000 * (Y - Z);
		return make_tuple(Lab_L, Lab_a, Lab_b);
	}

	static double R_revise(double x)
	{
		if (x > (6.0 / 29.0))
			x = pow(x, 3.0);
		else
			x = (x - 16.0 / 116.0) * 3.0 *pow(6.0 / 29.0, 2);
		return x;
	}

	static std::tuple<double, double, double> LabToRgb1(double l, double a, double b)
	{
		double  R, G, B = 0;
		//LAB->XYZ
		double var_Y = (l + 16) / 116.0;
		double var_X = var_Y + a / 500.0;
		double var_Z = var_Y - b / 200.0;

		var_Y = R_revise(var_Y);
		var_X = R_revise(var_X);
		var_Z = R_revise(var_Z);

		double X = 95.04 * var_X; 
		double Y = 100.000 * var_Y;
		double Z = 108.89 * var_Z;

		//XYZ->RGB
		var_X = X / 100.0;        //X from 0 to  95.047      (Observer = 2°, Illuminant = D65)
		var_Y = Y / 100.0;        //Y from 0 to 100.000
		var_Z = Z / 100.0;       //Z from 0 to 108.883

		double var_R = 3.24071 * var_X + (-1.53726) * var_Y + (-0.498571) * var_Z;
		double var_G = (-0.969258) * var_X + 1.87599 * var_Y + 0.0415557 * var_Z;
		double var_B = 0.0556352 * var_X + (-0.203996) * var_Y + 1.05707 * var_Z;

		if (var_R > 0.0031308)
			var_R = 1.055 * (pow(var_R, 1 / 2.4)) - 0.055;
		else
			var_R = 12.92 * var_R;

		if (var_G > 0.0031308)
			var_G = 1.055 * (pow(var_G, 1 / 2.4) - 0.055);
		else
			var_G = 12.92 * var_G;

		if (var_B > 0.0031308)
			var_B = 1.055 * (pow(var_B, 1 / 2.4) - 0.055);
		else
			var_B = 12.92 * var_B;

		R = (var_R)* 255.0;
		G = (var_G)* 255.0;
		B = (var_B)* 255.0;

		R = min(255.0, max(0.0, R));
		G = min(255.0, max(0.0, G));
		B = min(255.0, max(0.0, B));

		return make_tuple(R, G, B);
	}

	static std::vector<double> HSB2Lab(const double &h, const double & s, const double & b)
	{
		int rgb[3];
		double hsv[3] = {h,s,b};
		HSB2RGB(hsv, rgb);
		return RGB2Lab(rgb[0], rgb[1],rgb[2]);
	}

	static float getColorDH(const double& c1, const double& c2)
	{		
		if (c1 - c2 > 180)
			return 360 - c1 + c2;
		if (c2 - c1 > 180)
			return -(360 - c2 + c1);
		else
			return (c1 - c2);
	}

	static void HSB2RGB(const double hsv[3], int(&rgb)[3]) {
		
		
		int i;

		float RGB_min, RGB_max;
		RGB_max = hsv[2]*2.55f;
		RGB_min = RGB_max*(100 - hsv[1]) / 100.0f;

		i =hsv[0] / 60;
		int difs = (int)hsv[0] % 60; // factorial part of h

						   // RGB adjustment amount by hue 
		float RGB_Adj = (RGB_max - RGB_min)*difs / 60.0f;

		switch (i) {
		case 0:
			rgb[0] = RGB_max;
			rgb[1] = RGB_min + RGB_Adj;
			rgb[2] = RGB_min;
			break;
		case 1:
			rgb[0] = RGB_max - RGB_Adj;
			rgb[1] = RGB_max;
			rgb[2] = RGB_min;
			break;
		case 2:
			rgb[0] = RGB_min;
			rgb[1] = RGB_max;
			rgb[2] = RGB_min + RGB_Adj;
			break;
		case 3:
			rgb[0] = RGB_min;
			rgb[1] = RGB_max - RGB_Adj;
			rgb[2] = RGB_max;
			break;
		case 4:
			rgb[0] = RGB_min + RGB_Adj;
			rgb[1] = RGB_min;
			rgb[2] = RGB_max;
			break;
		default:		// case 5:
			rgb[0] = RGB_max;
			rgb[1] = RGB_min;
			rgb[2] = RGB_max - RGB_Adj;
			break;
		}

	
	}

	static std::vector<double> RGB2Lab(const double & R, const double & G, const double & B)
	{
		double Lab_L, Lab_a, Lab_b;
		Lab_L = Lab_a = Lab_b = 0;

		double X, Y, Z;
		double r = R / 255.000; // rgb range: 0 ~ 1
		double g = G / 255.000;
		double b = B / 255.000;
		// gamma 2.2
		if (r > 0.04045)
			r = pow((r + 0.055) / 1.055, 2.4);
		else
			r = r / 12.92;
		if (g > 0.04045)
			g = pow((g + 0.055) / 1.055, 2.4);
		else
			g = g / 12.92;
		if (b > 0.04045)
			b = pow((b + 0.055) / 1.055, 2.4);
		else
			b = b / 12.92;
		// sRGB
		X = r * 0.436052025 + g * 0.385081593 + b * 0.143087414;
		Y = r * 0.222491598 + g * 0.716886060 + b * 0.060621486;
		Z = r * 0.013929122 + g * 0.097097002 + b * 0.714185470;
		// XYZ range: 0~100
		X = X * 100.000;
		Y = Y * 100.000;
		Z = Z * 100.000;
		// Reference White Point
		double ref_X = 96.4221;
		double ref_Y = 100.000;
		double ref_Z = 82.5211;
		X = X / ref_X;
		Y = Y / ref_Y;
		Z = Z / ref_Z;
		// Lab
		if (X > 0.008856)
			X = pow(X, 1 / 3.000);
		else
			X = (7.787 * X) + (16 / 116.000);
		if (Y > 0.008856)
			Y = pow(Y, 1 / 3.000);
		else
			Y = (7.787 * Y) + (16 / 116.000);
		if (Z > 0.008856)
			Z = pow(Z, 1 / 3.000);
		else
			Z = (7.787 * Z) + (16 / 116.000);
		Lab_L = (116.000 * Y) - 16.000;
		Lab_a = 500.000 * (X - Y);
		Lab_b = 200.000 * (Y - Z);
		return{ Lab_L, Lab_a, Lab_b };
	}

	static double calColorDist(const vector<double> &l1, const vector<double> &l2)
	{
		return sqrt(pow(l1[0] - l2[0], 2) + pow(l1[1] - l2[1], 2) + pow(l1[2] - l2[2], 2));
	}

	static std::vector<double> RGB2HSV(double rgbR, double rgbG, double rgbV)
	{
		if (rgbR == rgbG && rgbG == rgbV)
		{
			rgbG++;
			rgbV--;
		}

		double hsbH, hsbS, hsbV;
		hsbH = hsbS = hsbV = 0.0;

		double _min = min(min(rgbR, rgbG), rgbV);
		double _max = max(max(rgbR, rgbG), rgbV);
		double del_Max = _max - _min;

		hsbV = _max / 255.0 * 100;
		hsbS = del_Max / _max * 100;

		if (_max == rgbR && rgbG >= rgbV)
		{
			hsbH = (rgbG - rgbV) * 60.0 / (_max - _min) + 0;
		}
		else if (_max == rgbR && rgbG < rgbV)
		{
			hsbH = (rgbG - rgbV) * 60.0 / (_max - _min) + 360;
		}
		else if (_max == rgbG)
		{
			hsbH = (rgbV - rgbR) * 60.0 / (_max - _min) + 120;
		}
		else if (_max == rgbV)
		{
			hsbH = (rgbR - rgbG) * 60.0 / (_max - _min) + 240;
		}
		return{ _min(360.0, _max(0.0, hsbH)) ,_min(100.0, _max(0.0, hsbS)) ,_min(100.0, _max(0.0, hsbV)) };
	}

	static std::vector<double> LabToRgb(double l, double a, double b)
	{
		double  R, G, B = 0;
		//LAB->XYZ
		double var_Y = (l + 16) / 116.0;
		double var_X = var_Y + a / 500.0;
		double var_Z = var_Y - b / 200.0;

		var_Y = R_revise(var_Y);
		var_X = R_revise(var_X);
		var_Z = R_revise(var_Z);

		double X = 95.04 * var_X;     //ps
		double Y = 100.000 * var_Y;
		double Z = 108.89 * var_Z;

		//XYZ->RGB
		var_X = X / 100.0;        //X from 0 to  95.047      (Observer = 2°, Illuminant = D65)
		var_Y = Y / 100.0;        //Y from 0 to 100.000
		var_Z = Z / 100.0;       //Z from 0 to 108.883

		double var_R = 3.24071 * var_X + (-1.53726) * var_Y + (-0.498571) * var_Z;
		double var_G = (-0.969258) * var_X + 1.87599 * var_Y + 0.0415557 * var_Z;
		double var_B = 0.0556352 * var_X + (-0.203996) * var_Y + 1.05707 * var_Z;

		if (var_R > 0.0031308)
			var_R = 1.055 * (pow(var_R, 1 / 2.4)) - 0.055;
		else
			var_R = 12.92 * var_R;

		if (var_G > 0.0031308)
			var_G = 1.055 * (pow(var_G, 1 / 2.4) - 0.055);
		else
			var_G = 12.92 * var_G;

		if (var_B > 0.0031308)
			var_B = 1.055 * (pow(var_B, 1 / 2.4) - 0.055);
		else
			var_B = 12.92 * var_B;

		R = (var_R)* 255.0;
		G = (var_G)* 255.0;
		B = (var_B)* 255.0;

		R = _min(255.0, _max(0.0, R));
		G = _min(255.0, _max(0.0, G));
		B = _min(255.0, _max(0.0, B));

		return{ R, G, B };
	}


private:

};



