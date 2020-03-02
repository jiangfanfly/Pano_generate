#pragma once
#include <opencv2\core\core.hpp>
#include <opencv2\highgui\highgui.hpp>
class myLadybugGeometry
{
public:
	myLadybugGeometry();
	~myLadybugGeometry();

public:
	//读入5个独立相机的内参保存在成员变量mvInnerParas中
	//内参读入顺序为：
	//CameraNum,k1,k2,k3,k4,p1,p2,c1,c2,lamda,fdistroted,x0distorted,y0distorted,frectified,x0rectified,y0rectified
	bool ReadInnerParaFromFile(std::string FilePath);

	bool ReadInvInnerParaFromFile(std::string FilePath);

	//读入5个独立相机的外参(外参的给出形式为：Rx,Ry,Rz,Tx,Ty,Tz,即三个角元素+三个线元素)
	//读入并转换为变换矩阵[R|t],保存在 mvUnitCamExParas中
	bool ReadUnitCameraExParaFromFile(std::string strFilePath);

	//将鱼眼图像的像素坐标转换为纠正后图像的像素坐标
	//输入第1参数为相机号，第2,3参数为输入的鱼眼像素坐标，第4,5参数为输出纠正后像素坐标
	bool LadybugRectifyImage(int CameraNum, double Pixalx, double Pixaly, 
		double* RectifiedPixalx, double* RectifiedPixaly);

	bool LadybugUnRectifyImage(int CameraNum, double Pixalx, double Pixaly, double* DistortedPixalx, double* DistortedPixaly);

	//将鱼眼图像的像素坐标投影到全景球面上,Radius是自定义的球半径(单位：m)
	bool LadybugProjectFishEyePtToSphere(int CameraNum, 
		double FishEyePixalx, double FishEyePixaly,double Radius,
		double* SphereX, double* SphereY, double* SphereZ);

	bool LadybugReprojectSpherePtToFishEyeImg(double SphereX, double SphereY, double SphereZ, double Radius,
		int* CameraNum, double* FishEyePixalx, double* FishEyePixaly);
	

	//将全景图像上的坐标投射到球面上;
	bool LadybugReprojectPanoPtToSphere(double SphereRadius, int PanoImgWith,int PanoImgHeight, double xPano, double yPano, 
		double* SphereX, double* SphereY, double* SphereZ);
	bool LadybugReprojectPanoPtToFishEyeImg(int PanoImgWith,int PanoImgHeight, double SphereRadius, double xPano, double yPano,
		int* CamID,double* FishEyePixalx, double* FishEyePixaly);

	//输入图像大小
	void InputImgSize(int width, int height);

	//输出纠正后图像的像素中心
	void GetRectifiedImgCenter(int CameraNum, double* x0, double *y0);
	
	//输出鱼眼图像的像素中心
	void GetFishEyeImgCenter(int CameraNum, double* x0, double* y0);

	//输出纠正后图像的焦距
	double GetRectifiedImgFocalLength(int CameraNum);

	//双线性内插函数
	cv::Vec3b BilinearInterpolation(cv::Mat Img, double x, double y);

private:
	double mImgHeight;
	double mImgWidth;
	std::vector<cv::Mat> mvInnerParas;
	std::vector<cv::Mat> mvInvInnerParas;
	std::vector<cv::Mat> mvUnitCamExParas;
	cv::Mat mCurrentInnerPara;
	cv::Mat mCurrentInvInnerPara;

	//由于Ladybug全景相机是横着安放的，其系统使用图像是横着的，我们习惯使用正常的图像(竖着的)
	//因此，操作时需要先从竖着的图像像素坐标转换为横着的像素坐标(逆时针旋转90度)
	//处理完成之后需要顺时针旋转90度转换回我们习惯的竖着的图像像素坐标
	cv::Point2d ImgPixCoor2LadybugPixCoor(cv::Point2d pt, int ImgWidth, int ImgHeight);
	cv::Point2d LadybugPixCoor2ImgPixCoor(cv::Point2d pt, int ImgWidth, int Height);

	/*通过给出的欧拉角Rx，Ry，Rz计算旋转矩阵，计算方式如下：	
	     |((cRz)(cRy)) ((cRz)(sRy)(sRx)-(sRz)(cRx)) ((cRz)(sRy)(cRx)+(sRz)(sRx)) Tx|
	 R = |((sRz)(cRy)) ((sRz)(sRy)(sRx)+(cRz)(cRx)) ((sRz)(sRy)(cRx)-(cRz)(sRx)) Ty|
	     |((-sRy))     ((cRy)(sRx))                 ((cRy)(cRx)))                Tz|
	     |0            0                            0                            1 |        */
	void CalculateRotationMatFromEulerAngle(double Rx, double Ry, double Rz, double* R);

	
	
};

