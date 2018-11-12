#ifndef HOUGHLINEFINDER_H
#define HOUGHLINEFINDER_H

#include <opencv2\highgui\highgui.hpp>
#include <opencv2\core\core.hpp>
#include <opencv2\imgproc\imgproc.hpp>
#include <vector>

#define M_PI 3.14159265358979323846

class HoughLineFinder
{
private:
	//输入图像
	cv::Mat img;
	//检测结果向量
	std::vector<cv::Vec4i> lines;
	//累加器分辨率参数
	double deltaRho;
	double deltaTheta;
	//结果最小投票数
	int minVote;
	//直线的最小长度
	double minLength;
	//直线上允许的最大空隙
	double maxGap;
public:
	//预设各参数值
	HoughLineFinder() :deltaRho(1), deltaTheta(M_PI / 180), minVote(10), minLength(0.0), maxGap(0.0) {}
	//设定分辨率参数
	void setAccResolution(double dRho, double dTheta)
	{
		deltaRho = dRho;
		deltaTheta = dTheta;
	}
	//设定最小投票数
	void setminVote(int minv)
	{
		minVote = minv;
	}
	//设定距离和间隔
	void setLengthAndGap(double length, double gap)
	{
		minLength = length;
		maxGap = gap;
	}
	//概率霍夫变换检测
	std::vector<cv::Vec4i> findLines(cv::Mat& binary)
	{
		lines.clear();
		cv::HoughLinesP(binary, lines, deltaRho, deltaTheta, minVote, minLength, maxGap);
		return lines;
	}
	//结果显示
	void drawDetectedLines(cv::Mat& image, cv::Scalar color = cv::Scalar(255, 255, 255))
	{
		std::vector<cv::Vec4i>::const_iterator it = lines.begin();
		while (it != lines.end())
		{
			cv::Point pt1((*it)[0], (*it)[1]);
			cv::Point pt2((*it)[2], (*it)[3]);
			cv::line(image, pt1, pt2, color);
			it++;
		}
	}
};
#endif