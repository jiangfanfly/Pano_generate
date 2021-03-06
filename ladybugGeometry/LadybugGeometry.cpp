// LadybugGeometry.cpp: 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "myLadybugGeometry.h"
#include <Windows.h>
#include <iostream>
#include <io.h>
#include <vector>


int CameraNum = 2;
const char strPath[] = "ImgFolderName";
std::string strFishEyeImgPath = "C:\\Users\\QZaneJ\\Desktop\\DC02.jpg";
std::string strFishEyeImgFolder = "E:\\02";
std::string strRectifiedImgPath = "C:\\Users\\QZaneJ\\Desktop\\RC02.jpg";
const char strInnerParaPath[] = "E:\\canshu\\jiejingInnPara.txt";
const char strExParaPath[] = "E:\\canshu\\jiejingExPara.txt";
const char strInvInnerParaPath[] = "E:\\canshu\\jiejingInvInnPara.txt";

const std::string strPathToSequences = "E:\\02";
std::string PanoImagePath = "E:\\02\\PanoImage2\\";
using namespace std;


void GetFileNames(std::string path, std::vector<std::string>& files)
{
	//文件句柄  
	__int64   hFile = 0;
	//文件信息  
	struct __finddata64_t fileinfo;
	std::string p;
	if ((hFile = _findfirst64(p.assign(path).append("\\*").c_str(), &fileinfo)) != -1)
	{
		std::cout << p << std::endl;
		do
		{
			//如果是目录,迭代之   
			//如果不是,加入列表  
			if ((fileinfo.attrib & _A_SUBDIR))
			{
				if (strcmp(fileinfo.name, ".") != 0 && strcmp(fileinfo.name, "..") != 0)
					GetFileNames(p.assign(path).append("\\").append(fileinfo.name), files);
			}
			else
			{
				files.push_back(p.assign(fileinfo.name));
			}
		} while (_findnext64(hFile, &fileinfo) == 0);
		_findclose(hFile);
	}
}
void LoadImagesForFishEye(const string &strPathToSequence, 
	vector<vector<string> > &vstrImageFishEye)
{
	
	vector<string> strPrefix00;
	vector<string> strPrefix01;
	vector<string> strPrefix02;
	vector<string> strPrefix03;
	vector<string> strPrefix04;
	vector<string> strPrefix05;


	GetFileNames(strPathToSequence + "\\Img0",strPrefix00);
	GetFileNames(strPathToSequence + "\\Img1",strPrefix01);
	GetFileNames(strPathToSequence + "\\Img2",strPrefix02);
	GetFileNames(strPathToSequence + "\\Img3",strPrefix03);
	GetFileNames(strPathToSequence + "\\Img4",strPrefix04);
	GetFileNames(strPathToSequence + "\\Img5",strPrefix05);


	const int nTimes = strPrefix00.size();
	cout << "Altogether: " << nTimes << endl;
	vstrImageFishEye.resize(nTimes);
	for (int i = 0; i<nTimes; i++)
	{
		vstrImageFishEye[i].resize(6);
	}

	for (int i = 0; i<nTimes; i++)
	{
		
		vstrImageFishEye[i][0] =  strPrefix00[i];
		vstrImageFishEye[i][1] =  strPrefix01[i];
		vstrImageFishEye[i][2] =  strPrefix02[i];
		vstrImageFishEye[i][3] =  strPrefix03[i];
		vstrImageFishEye[i][4] =  strPrefix04[i];
		vstrImageFishEye[i][5] =  strPrefix05[i];
	}

}
int main()
{
	/*std::vector<std::string> ImageFolder;
	ImageFolder.push_back("Img0");
	ImageFolder.push_back("Img1");
	ImageFolder.push_back("Img2");
	ImageFolder.push_back("Img3");
	ImageFolder.push_back("Img4");*/
	//int CameraNum = 2;
	myLadybugGeometry LBG;

	LBG.InputImgSize(1024,1224);
	bool isInnerParaLoaded = LBG.ReadInnerParaFromFile(strInnerParaPath);
	bool isInvInnerParaLoaded = LBG.ReadInvInnerParaFromFile(strInvInnerParaPath);
	bool isExParaLoaded = LBG.ReadUnitCameraExParaFromFile(strExParaPath);


	
	vector<vector<string> > strImageNames;
	LoadImagesForFishEye(strPathToSequences, strImageNames);
	

	/*std::string strFishEyeImgFolder_temp = strFishEyeImgFolder;
	FishEyeImg[0] = cv::imread(strFishEyeImgFolder.append("\\Img0\\").append("000001.jpg"),CV_LOAD_IMAGE_UNCHANGED);
	strFishEyeImgFolder = strFishEyeImgFolder_temp;
	FishEyeImg[1] = cv::imread(strFishEyeImgFolder.append("\\Img1\\").append("000001.jpg"), CV_LOAD_IMAGE_UNCHANGED);
	strFishEyeImgFolder = strFishEyeImgFolder_temp;
	FishEyeImg[2] = cv::imread(strFishEyeImgFolder.append("\\Img2\\").append("000001.jpg"), CV_LOAD_IMAGE_UNCHANGED);
	strFishEyeImgFolder = strFishEyeImgFolder_temp;
	FishEyeImg[3] = cv::imread(strFishEyeImgFolder.append("\\Img3\\").append("000001.jpg"), CV_LOAD_IMAGE_UNCHANGED);
	strFishEyeImgFolder = strFishEyeImgFolder_temp;
	FishEyeImg[4] = cv::imread(strFishEyeImgFolder.append("\\Img4\\").append("000001.jpg"), CV_LOAD_IMAGE_UNCHANGED);
	strFishEyeImgFolder = strFishEyeImgFolder_temp;
	FishEyeImg[5] = cv::imread(strFishEyeImgFolder.append("\\Img5\\").append("000001.jpg"), CV_LOAD_IMAGE_UNCHANGED);*/
	
	for (int i = 724, iend = strImageNames.size(); i != iend; i++)
	{
		std::cout << "Loading..." << std::endl;
		std::vector<cv::Mat> FishEyeImg(6);
		FishEyeImg[0] = cv::imread(strFishEyeImgFolder + "\\Img0\\" + strImageNames[i][0], CV_LOAD_IMAGE_UNCHANGED);
		FishEyeImg[1] = cv::imread(strFishEyeImgFolder + "\\Img1\\" + strImageNames[i][1], CV_LOAD_IMAGE_UNCHANGED);
		FishEyeImg[2] = cv::imread(strFishEyeImgFolder + "\\Img2\\" + strImageNames[i][2], CV_LOAD_IMAGE_UNCHANGED);
		FishEyeImg[3] = cv::imread(strFishEyeImgFolder + "\\Img3\\" + strImageNames[i][3], CV_LOAD_IMAGE_UNCHANGED);
		FishEyeImg[4] = cv::imread(strFishEyeImgFolder + "\\Img4\\" + strImageNames[i][4], CV_LOAD_IMAGE_UNCHANGED);
		FishEyeImg[5] = cv::imread(strFishEyeImgFolder + "\\Img5\\" + strImageNames[i][5], CV_LOAD_IMAGE_UNCHANGED);

		cv::Mat PanoImg = cv::Mat::zeros(2000, 4000, FishEyeImg[0].type());
		//cv::Mat PanoImg = cv::Mat::zeros(1000, 2000, FishEyeImg[0].type());
		std::cout << "Processing..." << std::endl;

		for (size_t i = 0; i < PanoImg.rows; i++)
		{
			for (size_t j = 0; j < PanoImg.cols; j++)
			{
				double xRectified, yRectified;
				int xDist = j;
				int yDist = i;
				int CamID;
				//LBG.LadybugReprojectPanoPtToFishEyeImg(2000, 1000, 20, xDist, yDist, &CamID, &xRectified, &yRectified);
				LBG.LadybugReprojectPanoPtToFishEyeImg(4000, 2000, 20, xDist, yDist, &CamID, &xRectified, &yRectified);
				if (xRectified >= 0 && xRectified < 1024 &&
					yRectified >= 0 && yRectified < 1224)
				{
					PanoImg.at<cv::Vec3b>(yDist, xDist) = LBG.BilinearInterpolation(FishEyeImg[CamID], xRectified, yRectified);
				}
				else if (i<1558)
				{
					if (xRectified < 0)
						xRectified = 0;
					if (yRectified < 0)
						yRectified = 0;
					if (xRectified >= 1024)
						xRectified = 1023;
					if (yRectified >= 1224)
						yRectified = 1223;

					PanoImg.at<cv::Vec3b>(yDist, xDist) = LBG.BilinearInterpolation(FishEyeImg[CamID], xRectified, yRectified);
				}
				else
				{
					PanoImg.at<cv::Vec3b>(yDist, xDist) = cv::Vec3b(0, 0, 0);
				}
			}

		}
		
		for (int i = 0; i < 6; i++)
		{
			FishEyeImg[i].release();
		}
		cv::imwrite(PanoImagePath+strImageNames[i][0], PanoImg);
		PanoImg.release();
		cout << "第 " << i+1 << " 张全景图像生成完成." << endl;
		cout << endl;
	}
	

	//double xFishEye1, yFishEye1, xFishEye2, yFishEye2;
	//double xRectified, yRectified;
	//double xDistorted, yDistorted;
	//LBG.GetFishEyeImgCenter(CameraNum, &xDistorted, &yDistorted);
	//double SphereX1, SphereY1, SphereZ1;
	//double SphereX2, SphereY2, SphereZ2;
	//bool flag = false;
	//xFishEye1 =xDistorted;
	//yFishEye1 = yDistorted;
	//xFishEye2 = xDistorted+1;
	//yFishEye2 = yDistorted+1;

	//double xFishEye3 = 0;
	//double yFishEye3 = 1616;

	//double x1 = 20;
	//double y1 = 0;
	//double z1 = 0;

	//std::cout << "SpherePoints: " << std::endl;
	//flag = LBG.LadybugProjectFishEyePtToSphere(1, xFishEye3, yFishEye3, 20, &SphereX1, &SphereY1, &SphereZ1);
	//std::cout << "mark:" << std::endl;
	//std::cout << SphereX1 << "  " << SphereY1 << "  " << SphereZ1 << std::endl;
	//LBG.LadybugProjectFishEyePtToSphere(1, xFishEye2, yFishEye2, 20, &SphereX2, &SphereY2, &SphereZ2);
	//std::cout << SphereX2 << "  " << SphereY2 << "  " << SphereZ2 << std::endl;
	//double DX = SphereX1 - SphereX2;
	//double DY = SphereY1 - SphereY2;
	//double DZ = SphereZ1 - SphereZ2;
	//std::cout << "BiasFromOnePixal:" << sqrt(DX*DX + DY * DY + DZ * DZ)<<std::endl;
	//LBG.LadybugProjectFishEyePtToSphere(1, xFishEye1, yFishEye2, 20, &SphereX1, &SphereY1, &SphereZ1);
	//std::cout << SphereX1 << "  " << SphereY1 << "  " << SphereZ1 << std::endl;
	//LBG.LadybugProjectFishEyePtToSphere(1, xFishEye2, yFishEye2, 20, &SphereX1, &SphereY1, &SphereZ1);
	//std::cout << SphereX1 << "  " << SphereY1 << "  " << SphereZ1 << std::endl;
	//std::cout << LBG.GetRectifiedImgFocalLength(4) << std::endl;

	//double rpx;
	//double rpy;
	//int Cam;

	//LBG.LadybugProjectFishEyePtToSphere(4, xFishEye3, yFishEye3, 20, &SphereX1, &SphereY1, &SphereZ1);
	//std::cout <<"rp: "<< SphereX1 << "  " << SphereY1 << "  " << SphereZ1 << std::endl;
	//LBG.LadybugReprojectSpherePtToFishEyeImg(SphereX1, SphereY1, SphereZ1, 20, &Cam, &rpx, &rpy);
	//std::cout<< "cam: " << Cam<<std::endl;
	//std::cout << "rp: " << rpx << " " << rpy << std::endl;

	//LBG.LadybugProjectFishEyePtToSphere(Cam, rpx, rpy, 20, &SphereX1, &SphereY1, &SphereZ1);
	//std::cout << "rp: " << SphereX1 << "  " << SphereY1 << "  " << SphereZ1 << std::endl;
	//getchar();
	//for (size_t n = 0; n <ImageFolder.size(); n++)
	//{
	//	CameraNum = n;
	//	std::string temp = strFishEyeImgPath;
	//	std::string temp2 = strFishEyeImgPath;
	//	temp.append("\\").append(ImageFolder[n]);
	//	temp2.append("\\").append("Rectified").append(ImageFolder[n]);
	//	std::vector<std::string> vImgNames;
	//	GetFileNames(temp, vImgNames);
	//	std::cout << "All together " << vImgNames.size() << " Images were found at " << strFishEyeImgPath << std::endl;
	//	std::cout << "Starting to LoadImage one by one and Rectify it...." << std::endl;
	//	std::cout << std::endl;

	//	myLadybugGeometry LBG;
	//	LBG.InputImgSize(1232, 1616);
	//	bool isInnerParaLoaded = LBG.ReadInnerParaFromFile(strInnerParaPath);
	//	bool isInvInnerParaLoaded = LBG.ReadInvInnerParaFromFile(strInvInnerParaPath);
	//	bool isExParaLoaded = LBG.ReadUnitCameraExParaFromFile(strExParaPath);
	//	for (size_t t = 0, tend = vImgNames.size(); t < tend; t++)
	//	{
	//		std::string p1 = temp;
	//		std::string p2 = temp2;
	//		std::string FishEyeImgPathName = p1.append("\\").append(vImgNames[t]);
	//		std::string RectifiedImgPathName = p2.append("\\").append(vImgNames[t]);
			//cv::Mat FishEyeImg = cv::imread(strFishEyeImgPath, 1);
			//cv::Mat RectifiedImg = cv::Mat::zeros(FishEyeImg.rows, FishEyeImg.cols, FishEyeImg.type());

			//for (size_t i = 0; i < RectifiedImg.rows; i++)
			//{
			//	for (size_t j = 0; j < RectifiedImg.cols; j++)
			//	{
			//		double xDistorted, yDistorted;
			//		int xRect = j;
			//		int yRect = i;
			//		int flag = LBG.LadybugUnRectifyImage(CameraNum, xRect, yRect, &xDistorted, &yDistorted);
			//		if (xDistorted >= 0 && xDistorted < RectifiedImg.cols &&
			//			yDistorted >= 0 && yDistorted < RectifiedImg.rows)
			//		{
			//			RectifiedImg.at<cv::Vec3b>(yRect, xRect) = FishEyeImg.at<cv::Vec3b>(yDistorted, xDistorted);
			//		}
			//	}

			//}
			//cv::imwrite(strRectifiedImgPath.c_str(), RectifiedImg);
			//RectifiedImg.release();
			//FishEyeImg.release();
	//		std::cout << "the " << vImgNames[t] << " has been rectified and saved !" << std::endl;
	//	}

	//}
	

	

	
	//cvNamedWindow("Rect");
	//cvShowImage("Rect", &IplImage(RectifiedImg));
	return 0;
}

