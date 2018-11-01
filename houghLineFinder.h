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
	//����ͼ��
	cv::Mat img;
	//���������
	std::vector<cv::Vec4i> lines;
	//�ۼ����ֱ��ʲ���
	double deltaRho;
	double deltaTheta;
	//�����СͶƱ��
	int minVote;
	//ֱ�ߵ���С����
	double minLength;
	//ֱ�������������϶
	double maxGap;
public:
	//Ԥ�������ֵ
	HoughLineFinder() :deltaRho(1), deltaTheta(M_PI / 180), minVote(10), minLength(0.0), maxGap(0.0) {}
	//�趨�ֱ��ʲ���
	void setAccResolution(double dRho, double dTheta)
	{
		deltaRho = dRho;
		deltaTheta = dTheta;
	}
	//�趨��СͶƱ��
	void setminVote(int minv)
	{
		minVote = minv;
	}
	//�趨����ͼ��
	void setLengthAndGap(double length, double gap)
	{
		minLength = length;
		maxGap = gap;
	}
	//���ʻ���任���
	std::vector<cv::Vec4i> findLines(cv::Mat& binary)
	{
		lines.clear();
		cv::HoughLinesP(binary, lines, deltaRho, deltaTheta, minVote, minLength, maxGap);
		return lines;
	}
	//�����ʾ
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