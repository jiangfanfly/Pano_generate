#include "stdafx.h"
#include "myLadybugGeometry.h"
#include <iostream>

myLadybugGeometry::myLadybugGeometry()
{
	
}

myLadybugGeometry::~myLadybugGeometry()
{
}

//读入5个独立相机的内参保存在成员变量mvInnerParas中
//内参读入顺序为：
//k1,k2,k3,k4,p1,p2,c1,c2,lamda,fdistroted,x0distorted,y0distorted,frectified,x0rectified,y0rectified
bool myLadybugGeometry::ReadInnerParaFromFile(std::string strFilePath)
{
	//相机编号
	int CameraNum;
	//畸变参数k1,k2,k3,k4,p1,p2,c1,c2
	double k1, k2, k3, k4;
	double p1, p2;
	double c1, c2;
	//尺度参数lamda
	double lamda;
	//鱼眼图像内定向参数(fdistorted, x0distorted, y0distorted)
	//矫正图像内定向参数(frectified，x0rectified, y0rectified)
	double fdistorted, frectified;
	double x0distorted, y0distorted;
	double x0rectified, y0rectified;
	
	//从文件中读入5个相机的内参
	FILE* fp = fopen(strFilePath.c_str(), "r");
	for (size_t i = 0; i < 6; i++)
	{
		mCurrentInnerPara = cv::Mat::zeros(cvSize(16, 1), CV_64FC1);
		fscanf(fp, "%d  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf\r\n",
			&CameraNum, &k1, &k2, &k3, &k4, &p1, &p2, &c1, &c2,
			&lamda, &fdistorted, &x0distorted, &y0distorted, &frectified, &x0rectified, &y0rectified);
		mCurrentInnerPara.at<double>(0, 0) = CameraNum;
		mCurrentInnerPara.at<double>(0, 1) = k1;
		mCurrentInnerPara.at<double>(0, 2) = k2;
		mCurrentInnerPara.at<double>(0, 3) = k3;
		mCurrentInnerPara.at<double>(0, 4) = k4;
		mCurrentInnerPara.at<double>(0, 5) = p1;
		mCurrentInnerPara.at<double>(0, 6) = p2;
		mCurrentInnerPara.at<double>(0, 7) = c1;
		mCurrentInnerPara.at<double>(0, 8) = c2;
		mCurrentInnerPara.at<double>(0, 9) = lamda;
		mCurrentInnerPara.at<double>(0, 10) = fdistorted;
		mCurrentInnerPara.at<double>(0, 11) = x0distorted;
		mCurrentInnerPara.at<double>(0, 12) = y0distorted;
		mCurrentInnerPara.at<double>(0, 13) = frectified;
		mCurrentInnerPara.at<double>(0, 14) = x0rectified;
		mCurrentInnerPara.at<double>(0, 15) = y0rectified;
		//分别压入mvInnerPara中
		mvInnerParas.push_back(mCurrentInnerPara);
		mCurrentInnerPara.release();
	}

	fclose(fp);
	if (mvInnerParas.size()==6)
	{
		return true;
	}
	else
	{
		std::cout << "内参读入有问题，请检查是否有5个相机的内参且每个相机有16个内参" << std::endl;
		return false;
	}
}

bool myLadybugGeometry::ReadInvInnerParaFromFile(std::string strFilePath)
{//相机编号
	int CameraNum;
	//畸变参数k1,k2,k3,k4,p1,p2,c1,c2
	double k1, k2, k3, k4;
	double p1, p2;
	double c1, c2;
	//尺度参数lamda
	double lamda;
	//鱼眼图像内定向参数(fdistorted, x0distorted, y0distorted)
	//矫正图像内定向参数(frectified，x0rectified, y0rectified)
	double fdistorted, frectified;
	double x0distorted, y0distorted;
	double x0rectified, y0rectified;

	//从文件中读入5个相机的内参
	FILE* fp = fopen(strFilePath.c_str(), "r");
	for (size_t i = 0; i < 6; i++)
	{
		mCurrentInvInnerPara = cv::Mat::zeros(cvSize(16, 1), CV_64FC1);
		fscanf(fp, "%d  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf\r\n",
			&CameraNum, &k1, &k2, &k3, &k4, &p1, &p2, &c1, &c2,
			&lamda, &fdistorted, &x0distorted, &y0distorted, &frectified, &x0rectified, &y0rectified);
		mCurrentInvInnerPara.at<double>(0, 0) = CameraNum;
		mCurrentInvInnerPara.at<double>(0, 1) = k1;
		mCurrentInvInnerPara.at<double>(0, 2) = k2;
		mCurrentInvInnerPara.at<double>(0, 3) = k3;
		mCurrentInvInnerPara.at<double>(0, 4) = k4;
		mCurrentInvInnerPara.at<double>(0, 5) = p1;
		mCurrentInvInnerPara.at<double>(0, 6) = p2;
		mCurrentInvInnerPara.at<double>(0, 7) = c1;
		mCurrentInvInnerPara.at<double>(0, 8) = c2;
		mCurrentInvInnerPara.at<double>(0, 9) = lamda;
		mCurrentInvInnerPara.at<double>(0, 10) = fdistorted;
		mCurrentInvInnerPara.at<double>(0, 11) = x0distorted;
		mCurrentInvInnerPara.at<double>(0, 12) = y0distorted;
		mCurrentInvInnerPara.at<double>(0, 13) = frectified;
		mCurrentInvInnerPara.at<double>(0, 14) = x0rectified;
		mCurrentInvInnerPara.at<double>(0, 15) = y0rectified;
		//分别压入mvInnerPara中
		mvInvInnerParas.push_back(mCurrentInvInnerPara);
		mCurrentInvInnerPara.release();
	}
	fclose(fp);
	if (mvInvInnerParas.size() == 6)
	{
		return true;
	}
	else
	{
		std::cout << "内参读入有问题，请检查是否有5个相机的内参且每个相机有16个内参" << std::endl;
		return false;
	}

}

//读入5个独立相机的外参 (外参的给出形式为：Rx,Ry,Rz,Tx,Ty,Tz,即三个角元素+三个线元素)
//读入并转换为变换矩阵[R|t],保存在 mvUnitCamExParas中
bool myLadybugGeometry::ReadUnitCameraExParaFromFile(std::string strFilePath)
{
	double Rx, Ry, Rz, Tx, Ty, Tz;
	double R[9];

	FILE* fp = fopen(strFilePath.c_str(), "r");
	for (size_t i = 0; i < 6; i++)
	{
		cv::Mat TransMat = cv::Mat::zeros(cvSize(4, 4), CV_64FC1);
		fscanf(fp, "%lf  %lf  %lf  %lf  %lf  %lf\r\n", &Rx, &Ry, &Rz, &Tx, &Ty, &Tz);
		//计算旋转矩阵
		CalculateRotationMatFromEulerAngle(Rx, Ry, Rz, R);
		TransMat.at<double>(0, 0) = R[0];
		TransMat.at<double>(0, 1) = R[1];
		TransMat.at<double>(0, 2) = R[2];
		TransMat.at<double>(0, 3) = Tx;
		TransMat.at<double>(1, 0) = R[3];
		TransMat.at<double>(1, 1) = R[4];
		TransMat.at<double>(1, 2) = R[5];
		TransMat.at<double>(1, 3) = Ty;
		TransMat.at<double>(2, 0) = R[6];
		TransMat.at<double>(2, 1) = R[7];
		TransMat.at<double>(2, 2) = R[8];
		TransMat.at<double>(2, 3) = Tz;
		TransMat.at<double>(3, 0) = 0.0;
		TransMat.at<double>(3, 1) = 0.0;
		TransMat.at<double>(3, 2) = 0.0;
		TransMat.at<double>(3, 3) = 1.0;
		//压入成员变量mvUnitCamExParas变量中
		mvUnitCamExParas.push_back(TransMat.clone());
		TransMat.release();
	}
	if (mvUnitCamExParas.size()==6)
	{
		return true;

	}
	else
	{
		std::cout << "请检查外参文件格式(5个相机，每个相机有Rx,Ry,Rz,Tx,Ty,Tz)" << std::endl;
		return false;
	}
	
}

//将鱼眼图像的像素坐标转换为纠正后图像的像素坐标
//输入第1参数为相机号，第2,3参数为输入的鱼眼像素坐标，第4,5参数为输出纠正后像素坐标
bool myLadybugGeometry::LadybugRectifyImage(int CameraNum,
	double Pixalx, double Pixaly, 
	double* RectifiedPixalx, double* RectifiedPixaly)
{
	//畸变参数k1,k2,k3,k4,p1,p2,c1,c2
	double k1, k2, k3, k4;
	double p1, p2;
	double c1, c2;
	//尺度参数lamda
	double lamda;
	//鱼眼图像内定向参数(fdistorted, x0distorted, y0distorted)
	//矫正图像内定向参数(frectified，x0rectified, y0rectified)
	double fdistorted, frectified;
	double x0distorted, y0distorted;
	double x0rectified, y0rectified;

	mCurrentInnerPara = mvInnerParas[CameraNum];
	if (CameraNum != int(mCurrentInnerPara.at<double>(0,0)))
	{
		for (size_t i = 0; i < mvInnerParas.size(); i++)
		{
			if (CameraNum = int(mvInnerParas[i].at<double>(0, 0)))
			{
				mCurrentInnerPara = mvInnerParas[i].clone();
			}
		}
	}

	//赋值得到每个内参参数
	k1 = mCurrentInnerPara.at<double>(0, 1);
	k2 = mCurrentInnerPara.at<double>(0, 2);
	k3 = mCurrentInnerPara.at<double>(0, 3);
	k4 = mCurrentInnerPara.at<double>(0, 4);
	p1 = mCurrentInnerPara.at<double>(0, 5);
	p2 = mCurrentInnerPara.at<double>(0, 6);
	c1 = mCurrentInnerPara.at<double>(0, 7);
	c2 = mCurrentInnerPara.at<double>(0, 8);
	lamda = mCurrentInnerPara.at<double>(0, 9);
	fdistorted = mCurrentInnerPara.at<double>(0, 10);
	x0distorted = mCurrentInnerPara.at<double>(0, 11);
	y0distorted = mCurrentInnerPara.at<double>(0, 12);
	frectified = mCurrentInnerPara.at<double>(0, 13);
	x0rectified = mCurrentInnerPara.at<double>(0, 14);
	y0rectified = mCurrentInnerPara.at<double>(0, 15);

	//开始纠正
	//首先转换为LadyBug坐标系(相当于图像坐标系逆时针旋转90度)下的平面坐标，此处使用的是鱼眼相机的x0,y0;
	cv::Point2d LadybugPixalCoor = ImgPixCoor2LadybugPixCoor(cv::Point2d(Pixalx, Pixaly), mImgWidth, mImgHeight);
	double xt = LadybugPixalCoor.x- x0distorted;
	double yt = LadybugPixalCoor.y - y0distorted;
	//使用等距投影模型进行投影将鱼眼图像投影到纠正后图像
	double sqrtr = sqrt(xt*xt + yt * yt);
	double Rectifiedx = ((fdistorted*xt*tan(sqrtr / fdistorted)) / sqrtr);
	double Rectifiedy = ((fdistorted*yt*tan(sqrtr / fdistorted)) / sqrtr);
	//按照缩放比例缩放至鱼眼相机同等大小的尺度上
	double x_ = lamda * Rectifiedx;
	double y_ = lamda * Rectifiedy;
	//使用图像畸变模型，计算畸变改正量，使用的畸变参数有k1,k2,k3,k4,p1,p2,c1,c2
	double r2 = x_ * x_ + y_ * y_;
	double dx = x_ * (k1*r2 + k2 * r2*r2 + k3 * r2*r2*r2 + k4 * r2*r2*r2*r2) + 2 * p1*x_*y_ + p2 * (r2 + 2 * x_*x_) + c1 * x_ + c2 * y_;
	double dy = y_ * (k1*r2 + k2 * r2*r2 + k3 * r2*r2*r2 + k4 * r2*r2*r2*r2) + p1 * (r2 + 2 * y_*y_) + 2 * p2*x_*y_ + c1 * y_ + c2 * x_;
	//进行畸变改正并转换为像素坐标,注意此处要使用纠正后影像的x0,y0;
	double xRectTemp= (x_ - dx) + x0rectified;
	double yRectTemp= (y_ - dy) + y0rectified;
	//从Ladybug的像素坐标转换为正常图像的像素坐标(相当于顺时针旋转90度)
	cv::Point2d RectifiedPixalCoor = LadybugPixCoor2ImgPixCoor(cv::Point2d(xRectTemp, yRectTemp), mImgWidth, mImgHeight);
	*RectifiedPixalx = RectifiedPixalCoor.x;
	*RectifiedPixaly = RectifiedPixalCoor.y;

	return true;
}

bool myLadybugGeometry::LadybugUnRectifyImage(int CameraNum, double Pixalx, double Pixaly, double * DistortedPixalx, double * DistortedPixaly)
{
	//畸变参数k1,k2,k3,k4,p1,p2,c1,c2
	double k1, k2, k3, k4;
	double p1, p2;
	double c1, c2;
	//尺度参数lamda
	double lamda;
	//鱼眼图像内定向参数(fdistorted, x0distorted, y0distorted)
	//矫正图像内定向参数(frectified，x0rectified, y0rectified)
	double fdistorted, frectified;
	double x0distorted, y0distorted;
	double x0rectified, y0rectified;

	mCurrentInvInnerPara = mvInvInnerParas[CameraNum];
	if (CameraNum != int(mCurrentInvInnerPara.at<double>(0, 0)))
	{
		for (size_t i = 0; i < mvInvInnerParas.size(); i++)
		{
			if (CameraNum = int(mvInvInnerParas[i].at<double>(0, 0)))
			{
				mCurrentInnerPara = mvInvInnerParas[i].clone();
			}
		}
	}

	//赋值得到每个内参参数
	k1 = mCurrentInvInnerPara.at<double>(0, 1);
	k2 = mCurrentInvInnerPara.at<double>(0, 2);
	k3 = mCurrentInvInnerPara.at<double>(0, 3);
	k4 = mCurrentInvInnerPara.at<double>(0, 4);
	p1 = mCurrentInvInnerPara.at<double>(0, 5);
	p2 = mCurrentInvInnerPara.at<double>(0, 6);
	c1 = mCurrentInvInnerPara.at<double>(0, 7);
	c2 = mCurrentInvInnerPara.at<double>(0, 8);
	lamda = mCurrentInvInnerPara.at<double>(0, 9);
	fdistorted = mCurrentInvInnerPara.at<double>(0, 10);
	x0distorted = mCurrentInvInnerPara.at<double>(0, 11);
	y0distorted = mCurrentInvInnerPara.at<double>(0, 12);
	frectified = mCurrentInvInnerPara.at<double>(0, 13);
	x0rectified = mCurrentInvInnerPara.at<double>(0, 14);
	y0rectified = mCurrentInvInnerPara.at<double>(0, 15);

	cv::Point2d LadybugPixalCoor = ImgPixCoor2LadybugPixCoor(cv::Point2d(Pixalx, Pixaly), mImgWidth, mImgHeight);
	double x_ = LadybugPixalCoor.x - x0rectified;
	double y_ = LadybugPixalCoor.y - y0rectified;
	double r2 = x_ * x_ + y_ * y_;
	double dx = x_ * (k1*r2 + k2 * r2*r2 + k3 * r2*r2*r2 + k4 * r2*r2*r2*r2) + 2 * p1*x_*y_ + p2 * (r2 + 2 * x_*x_) + c1 * x_ + c2 * y_;
	double dy = y_ * (k1*r2 + k2 * r2*r2 + k3 * r2*r2*r2 + k4 * r2*r2*r2*r2) + p1 * (r2 + 2 * y_*y_) + 2 * p2*x_*y_ + c1 * y_ + c2 * x_;
	double xt = x_ + dx;
	double yt = y_ + dy;
	double sqrtr = sqrt(xt*xt + yt * yt);
	double Distortedx = ((frectified*xt*atan(sqrtr / frectified)) / sqrtr);
	double Distortedy = ((frectified*yt*atan(sqrtr / frectified)) / sqrtr);
	double xd = lamda*Distortedx + x0distorted;
	double yd = lamda*Distortedy + y0distorted;

	cv::Point2d DistortedPixalCoor = LadybugPixCoor2ImgPixCoor(cv::Point2d(xd, yd), mImgWidth, mImgHeight);
	*DistortedPixalx = DistortedPixalCoor.x;
	*DistortedPixaly = DistortedPixalCoor.y;

	return false;
}

//将鱼眼图像的像素坐标投影到全景球面上，Radius是自定义的球半径(单位：m)
bool myLadybugGeometry::LadybugProjectFishEyePtToSphere(int CameraNum, 
	double FishEyePixalx, double FishEyePixaly, double Radius,
	double* SphereX, double* SphereY, double* SphereZ)
{
	*SphereX = 0;
	*SphereY = 0;
	*SphereZ = 0;

	double RectifiedPixalx, RectifiedPixaly;
	cv::Point2d RectifiedPixalCoor;
	//从鱼眼图像像素坐标转换到纠正像素坐标
	LadybugRectifyImage(CameraNum, FishEyePixalx, FishEyePixaly, &RectifiedPixalx, &RectifiedPixaly);
	RectifiedPixalCoor = ImgPixCoor2LadybugPixCoor(cv::Point2d(RectifiedPixalx, RectifiedPixaly), mImgWidth, mImgHeight);
	RectifiedPixalx = RectifiedPixalCoor.x;
	RectifiedPixaly = RectifiedPixalCoor.y;

	//获取当前相机的内外方位元素
	cv::Mat CurrentTransMat = mvUnitCamExParas[CameraNum].clone();
	double FocalLength = GetRectifiedImgFocalLength(CameraNum);
	double x0rectified, y0rectified;
    GetRectifiedImgCenter(CameraNum,&x0rectified,&y0rectified);
	//得到的x0，y0使我们习惯的竖着的图像的像素中心，在Ladybug框架下计算需要换算成Ladybug框架下的像素中心
	cv::Point2d ImgCenter = ImgPixCoor2LadybugPixCoor(cv::Point2d(x0rectified, y0rectified), mImgWidth, mImgHeight);
	x0rectified = ImgCenter.x;
	y0rectified = ImgCenter.y;
	

	//组成此点在此相机坐标系下的坐标
	cv::Mat RectifiedImgPtCoor = cv::Mat::zeros(cvSize(1, 3), CV_64FC1);
	RectifiedImgPtCoor.at<double>(0, 0) = RectifiedPixalx - x0rectified;
	RectifiedImgPtCoor.at<double>(1, 0) = RectifiedPixaly - y0rectified;
	RectifiedImgPtCoor.at<double>(2, 0) = FocalLength;

	//获取当前相机的旋转矩阵与平移量
	cv::Mat RatationMat = CurrentTransMat.rowRange(0, 3).colRange(0, 3).clone();
	double Tx = CurrentTransMat.at<double>(0, 3);
	double Ty = CurrentTransMat.at<double>(1, 3);
	double Tz = CurrentTransMat.at<double>(2, 3);

	//此点乘以旋转矩阵以计算出其在Ladybug坐标系下的坐标，即在球面坐标系下的光线向量
	cv::Mat RayCoor = RatationMat * RectifiedImgPtCoor;
	double CX = RayCoor.at<double>(0, 0);
	double CY = RayCoor.at<double>(1, 0);
	double CZ = RayCoor.at<double>(2, 0);

	//计算向量长度，即光线与球面交会，球面半径为输入的Radius
	double a = CX * CX + CY * CY + CZ * CZ;
	double b = 2 * (CX * Tx + CY * Ty + CZ * Tz);
	double c = Tx * Tx + Ty * Ty + Tz * Tz - Radius * Radius;
	double b4ac = sqrt(b*b - 4 * a*c);
	double k = (-b + b4ac) / (2 * a);
	
	//其在球面下的坐标等于相机投影中心的位置向量+光线向量
	*SphereX = k * CX + Tx;
	*SphereY = k * CY + Ty;
	*SphereZ = k * CZ + Tz;

	if (SphereX != 0 || SphereY != 0 || SphereY != 0)
	{
		return true;
	}
	else
	{
		std::cout << "鱼眼图像坐标转球面坐标时出现问题，具体什么问题自己检查吧~嘻嘻" << std::endl;
		return false;
	}	
}

bool myLadybugGeometry::LadybugReprojectSpherePtToFishEyeImg(double SphereX, double SphereY, double SphereZ, double Radius, int * CameraNum, double * FishEyePixalx, double * FishEyePixaly)
{
	int Cam;
	cv::Mat T;
	cv::Mat RaySphere(3, 1, CV_64F) ;
	RaySphere.at<double>(0, 0) = SphereX;
	RaySphere.at<double>(1, 0) = SphereY;
	RaySphere.at<double>(2, 0) = SphereZ;
	double minCosAngle = 0.0;
	for (int i = 0; i < 6; i++)
	{
		cv::Mat CamPose=mvUnitCamExParas[i].col(3).rowRange(0,3);
		double cosAngle = double(RaySphere.dot(CamPose)) / ((cv::norm(RaySphere))*(cv::norm(CamPose)));
		if (cosAngle>0 && cosAngle <=1 && cosAngle>minCosAngle)
		{
			minCosAngle = cosAngle;
			T = mvUnitCamExParas[i].clone();
			Cam = i;
		}
		CamPose.release();
	}
	if (T.empty())
		return false;

	cv::Mat RotationMat = T.rowRange(0, 3).colRange(0, 3).clone();
	cv::Mat CamPos = T.col(3).rowRange(0, 3);
	cv::Mat RayCoor = RaySphere - CamPos;
	cv::Mat RectifiedImgPtCoor = RotationMat.inv()*RayCoor;
	double x = RectifiedImgPtCoor.at<double>(0, 0);
	double y = RectifiedImgPtCoor.at<double>(1, 0);
	double z = RectifiedImgPtCoor.at<double>(2, 0);

	double FocalLength = GetRectifiedImgFocalLength(Cam);
	double x0rectified, y0rectified;
	GetRectifiedImgCenter(Cam, &x0rectified, &y0rectified);

	cv::Point2d ImgCenter = ImgPixCoor2LadybugPixCoor(cv::Point2d(x0rectified, y0rectified), mImgWidth, mImgHeight);
	x0rectified = ImgCenter.x;
	y0rectified = ImgCenter.y;

	double RectPixalx = FocalLength * (x / z) + x0rectified;
	double RectPixaly = FocalLength * (y / z) + y0rectified;

	cv::Point2d RectPixalCoor = LadybugPixCoor2ImgPixCoor(cv::Point2d(RectPixalx, RectPixaly), mImgWidth, mImgHeight);
	RectPixalx = RectPixalCoor.x;
	RectPixaly = RectPixalCoor.y;
	double DistortedPixalx, DistortedPixaly;
	LadybugUnRectifyImage(Cam, RectPixalx, RectPixaly, &DistortedPixalx, &DistortedPixaly);
	*CameraNum = Cam;
	*FishEyePixalx = DistortedPixalx;
	*FishEyePixaly = DistortedPixaly;

	return true;
}

bool myLadybugGeometry::LadybugReprojectPanoPtToSphere(double SphereRadius, int PanoImgWith, int PanoImgHeight, 
	double xPano, double yPano, 
	double * SphereX, double * SphereY, double * SphereZ)
{
	double phh = CV_PI * (PanoImgWith - (2 * xPano)) / PanoImgWith;
	double phv = CV_PI * (PanoImgHeight - 2 * yPano) / (2 * PanoImgHeight);
	
	*SphereX = SphereRadius * cos(phv)*cos(phh);
	*SphereY = SphereRadius * cos(phv)*sin(phh);
	*SphereZ = SphereRadius * sin(phv);

	return true;
}

bool myLadybugGeometry::LadybugReprojectPanoPtToFishEyeImg(int PanoImgWidth,int PanoImgHeight,double SphereRadius, double xPano, double yPano,
	int* CamID,double* FishEyePixalx, double * FishEyePixaly)
{
	double tempSphereX, tempSphereY, tempSphereZ;
	double tempFishEyePixalx, tempFishEyePixaly;
	int tempCamID;
	LadybugReprojectPanoPtToSphere(SphereRadius,PanoImgWidth,PanoImgHeight,
		xPano,yPano,&tempSphereX,&tempSphereY,&tempSphereZ);
	LadybugReprojectSpherePtToFishEyeImg(tempSphereX, tempSphereY, tempSphereZ, SphereRadius, 
		&tempCamID, &tempFishEyePixalx, &tempFishEyePixaly);
	*CamID = tempCamID;
	*FishEyePixalx = tempFishEyePixalx;
	*FishEyePixaly = tempFishEyePixaly;
	return false;
}



//输入图像大小
void myLadybugGeometry::InputImgSize(int width, int height)
{
	mImgHeight = height;
	mImgWidth = width;
}

void myLadybugGeometry::GetRectifiedImgCenter(int CameraNum, double * x0, double * y0)
{
	mCurrentInnerPara = mvInnerParas[CameraNum];
	cv::Point2d ImgCenter = LadybugPixCoor2ImgPixCoor(cv::Point2d
	(mCurrentInnerPara.at<double>(0, 14), mCurrentInnerPara.at<double>(0, 15)),
		mImgWidth, mImgHeight);
	 *x0 = ImgCenter.x;
	 *y0 = ImgCenter.y;
}

void myLadybugGeometry::GetFishEyeImgCenter(int CameraNum, double * x0, double * y0)
{
	mCurrentInnerPara = mvInnerParas[CameraNum];
	cv::Point2d ImgCenter = LadybugPixCoor2ImgPixCoor(cv::Point2d
	(mCurrentInnerPara.at<double>(0, 11), mCurrentInnerPara.at<double>(0, 12)),
		mImgWidth, mImgHeight);
	*x0 = ImgCenter.x;
	*y0 = ImgCenter.y;
}

double myLadybugGeometry::GetRectifiedImgFocalLength(int CameraNum)
{
	mCurrentInnerPara = mvInnerParas[CameraNum];
	return mCurrentInnerPara.at<double>(0, 13);
}

cv::Point2d myLadybugGeometry::ImgPixCoor2LadybugPixCoor(cv::Point2d pt, int ImgWidth, int ImgHeight)
{
	double xImg, yImg, xLB, yLB;
	xImg = pt.x;
	yImg = pt.y;
	xLB = yImg;
	yLB = ImgWidth - xImg - 1;
	return cv::Point2d(xLB, yLB);

}

cv::Point2d myLadybugGeometry::LadybugPixCoor2ImgPixCoor(cv::Point2d pt, int ImgWidth, int Height)
{
	double xImg, yImg, xLB, yLB;
	xLB = pt.x;
	yLB = pt.y;
	xImg = ImgWidth - yLB - 1;
	yImg = xLB;
	return cv::Point2d(xImg, yImg);

}

/*通过给出的欧拉角Rx，Ry，Rz计算旋转矩阵，计算方式如下：
	|((cRz)(cRy)) ((cRz)(sRy)(sRx)-(sRz)(cRx)) ((cRz)(sRy)(cRx)+(sRz)(sRx)) Tx|
R = |((sRz)(cRy)) ((sRz)(sRy)(sRx)+(cRz)(cRx)) ((sRz)(sRy)(cRx)-(cRz)(sRx)) Ty|
	|((-sRy))     ((cRy)(sRx))                 ((cRy)(cRx)))                Tz|
	|0            0                            0                            1 |        */
void myLadybugGeometry::CalculateRotationMatFromEulerAngle(double Rx, double Ry, double Rz, double * R)
{
	double cRx, cRy, cRz, sRx, sRy, sRz;
	cRx = cos(Rx); cRy = cos(Ry); cRz = cos(Rz);
	sRx = sin(Rx); sRy = sin(Ry); sRz = sin(Rz);
	R[0] = cRz * cRy;
	R[1] = cRz * sRy * sRx - sRz * cRx;
	R[2] = cRz * sRy * cRx + sRz * sRx;
	R[3] = sRz * cRy;
	R[4] = sRz * sRy* sRx + cRz * cRx;
	R[5] = sRz * sRy* cRx - cRz * sRx;
	R[6] = -sRy;
	R[7] = cRy * sRx;
	R[8] = cRy * cRx;
}

cv::Vec3b myLadybugGeometry::BilinearInterpolation(cv::Mat img, double x, double y)
{
	int imgHeight = img.rows;
	int imgWidth = img.cols;
	int x11, x12, x21, x22, y11, y12, y21, y22;
	double dx, dy;
	cv::Vec3b I11, I12, I21, I22, I;
	x11 = int(x);
	y11 = int(y);
	x12 = x11 + 1;
	y12 = y11;
	x21 = x11;
	y21 = y11 + 1;
	x22 = x11 + 1;
	y22 = y11 + 1;
	dx = x - x11;
	dy = y - y11;
	if (x >= 0 && y >= 0 && x<imgWidth - 1 && y<imgHeight - 1)
	{
		I11 = img.at<cv::Vec3b>(y11, x11);
		I12 = img.at<cv::Vec3b>(y12, x12);
		I21 = img.at<cv::Vec3b>(y21, x21);
		I22 = img.at<cv::Vec3b>(y22, x22);
		I = (1 - dy)*(1 - dx)*I11 + (1 - dy)*dx*I12 + dy * (1 - dx)*I21 + dx * dy*I22;
	}
	else if (int(x) == imgWidth - 1 || int(y) == imgHeight - 1)
	{
		I = img.at<cv::Vec3b>(int(y), int(x));
	}
	else
	{
		I = cv::Vec3b(0,0,0);
	}
	return I;
}