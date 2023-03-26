#ifndef _PictureOP_H
#define _PictureOP_H
#pragma once

#include <vector>
#include <string>
#include<io.h>
#include<direct.h>

using namespace std;
class PictureOP
{
private:
	


public:
	static string getImageName(const std::string &picturepath)
	{
		char szDrive[_MAX_DRIVE];   
		char szDir[_MAX_DIR];       
		char szFname[_MAX_FNAME];   
		char szExt[_MAX_EXT];       
		_splitpath_s(picturepath.c_str(), szDrive, szDir, szFname, szExt);

		return szFname;
	}

	static string getMapName(const std::string &mappath)
	{
		
		char szDrive[_MAX_DRIVE];   
		char szDir[_MAX_DIR];       
		char szFname[_MAX_FNAME];   
		char szExt[_MAX_EXT];       
		_splitpath_s(mappath.c_str(), szDrive, szDir, szFname, szExt);

		string name = szFname;
		string xx = "_mapfile";
		return name.substr(0, name.size() - xx.size());
	}

	static bool exists(const std::string& name) 
	{
		return (_access(name.c_str(), 0) != -1);
	}

	static void getBackgroundPicture(const std::string &picturepath,std::string &backpath)
	{
		char szDrive[_MAX_DRIVE];   
		char szDir[_MAX_DIR];      
		char szFname[_MAX_FNAME];   
		char szExt[_MAX_EXT];      
		_splitpath_s(picturepath.c_str(), szDrive, szDir, szFname, szExt);

		string ImagefoldName;
		ImagefoldName = "extractedFeature\\background\\";

		backpath =  ImagefoldName + std::string(szFname) + "_back" + std::string(szExt);
	}

	

	static void getRelationshipPath(const std::string &mappath, std::string &Relationshippath)
	{
		char szDrive[_MAX_DRIVE];   
		char szDir[_MAX_DIR];       
		char szFname[_MAX_FNAME];   
		char szExt[_MAX_EXT];       
		_splitpath_s(mappath.c_str(), szDrive, szDir, szFname, szExt);

		
		Relationshippath = std::string(szDrive) + std::string(szDir) + std::string(szFname) + "_relationship" + std::string(szExt);
		if (_access(Relationshippath.c_str(), 0) == 0)
		{

		}
		else
		{
			//Relationshippath = "";
		}
	}
};
#endif
