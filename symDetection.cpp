//对一幅图像进行二值化

#include <iostream>
#include <opencv2\highgui\highgui.hpp>
#include <opencv2\core\core.hpp>
#include <opencv2\imgproc\imgproc.hpp>
#include <math.h>
#include <fstream>
#include "Harris.h"

#define M_PI 3.14159265358979323846
#define INFINITE_NUMBER 10000



using namespace std;
using namespace cv;

const int w = 500;
int levels = 0;
vector<vector<Point> > contours;
vector<Vec4i> hierarchy;


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

bool one_lines(vector<Point> line1, vector<Point> line2)
{
	Vec4f l1, l2;
	fitLine(line1, l1, CV_DIST_L2, 0, 0.01, 0.01);
	fitLine(line2, l2, CV_DIST_L2, 0, 0.01, 0.01);
	double k1 = l1[1] / l1[0];
	double k2 = l2[1] / l2[0];
	double b1 = l1[3] - k1 * l1[2];
	double b2 = l2[3] - k2 * l2[2];
	double dis = abs(b1 - b2) / sqrt(k1 * k1 + 1);
	if (dis < 0.00001)
		return true;
	else
		return false;
}

bool distance_equal(vector<Point> line1, vector<Point> line2)
{
	double len1 = sqrt(abs((line1[0].x - line1[1].x) * (line1[0].x - line1[1].x)) + abs((line1[0].y - line1[1].y) * (line1[0].y - line1[1].y)));
	double len2 = sqrt(abs((line2[0].x - line2[1].x) * (line2[0].x - line2[1].x)) + abs((line2[0].y - line2[1].y) * (line2[0].y - line2[1].y)));
	if (abs(len1 - len2) < 5)
		return true;
	else
		return false;
}

int min_x(Point p1, Point p2)
{
	if (p1.x < p2.x)
		return p1.x;
	else
		return p2.x;
}
int min_y(Point p1, Point p2)
{
	if (p1.y < p2.y)
		return p1.y;
	else		  
		return p2.y;
}
int max_x(Point p1, Point p2)
{
	if (p1.x > p2.x)
		return p1.x;
	else
		return p2.x;
}
int max_y(Point p1, Point p2)
{
	if (p1.y > p2.y)
		return p1.y;
	else
		return p2.y;
}

void draw_symaxis(Mat img, Point l1_p1, Point l1_p2, Point l2_p1, Point l2_p2)
{
	Point result_p1, result_p2;
	result_p1.x = (min_x(l1_p1, l1_p2) + min_x(l2_p1, l2_p2)) / 2;
	result_p1.y = (min_y(l1_p1, l1_p2) + min_y(l2_p1, l2_p2)) / 2;
	result_p2.x = (max_x(l1_p1, l1_p2) + max_x(l2_p1, l2_p2)) / 2;
	result_p2.y = (max_y(l1_p1, l1_p2) + max_y(l2_p1, l2_p2)) / 2;
	line(img, result_p1, result_p2, Scalar::all(125));
}

//寻找最下端点和最上端点
vector<Point> get_boud_point( vector<Point> line1, vector<Point> line2)
{
	vector<Point> result;
	return result;
}

//检测对称轴
void detectsymaxis(Mat img, vector<Point> corners, int r, int c)
{
	int size = corners.size();
	vector<Point>::iterator it;
	for (it = corners.begin(); it != corners.end(); it++)
	{
		it = corners.begin();
		if (it->x <= 2 || it->x >= c - 5 || it->y <= 2 || it->y >= r - 10)
		{
			corners.erase(it);
			it = corners.begin();
		}
		if (corners.size() == size - 4)
			break;
	}
	cout << "phase 7" << endl;
	//ofstream fout("1.txt");
	vector<double> theta;
	for (int i = 0; i < corners.size() - 1; i++)
	{
		double dx = corners[i].x - corners[i + 1].x;
		double dy = corners[i].y - corners[i + 1].y;
		if (dx == 0)
		{
			theta.push_back(INFINITE_NUMBER);
			//fout << INFINITE_NUMBER << endl;
		}
		else if (dy == 0)
		{
			theta.push_back(0);
			//fout << 0 << endl;
		}
		else
		{
			theta.push_back(dy / dx);
			//fout << abs(dy / dx) << endl;
		}
	}
	//fout.close();
	cout << "phase 8" << endl;
	int num = 0;
	for (int j = 0; j < theta.size() - 1; j++)
	{
		for (int i = 0; i < theta.size() - 1; i++)
		{
			if (abs(theta[j] - theta[i + 1]) < 0.8)
			{
				double b1 = corners[j].y - theta[j] * corners[j].x;
				double b2 = corners[i + 1].y - theta[i + 1] * corners[i + 1].x;
				double dis = abs(b1 - b2) / sqrt(theta[j] * theta[j] + 1);
				
				if (dis < 0.0001)
					continue;
				else
				{
					if (distance_equal)
					{
						draw_symaxis(img, corners[i], corners[i + 1], corners[j], corners[j + 1]);
						num++;
					}
					else
						continue;
				}

			}
			else
				continue;
		}
	}
	cout << "final phase" << endl;
	cout << "Number of symmetry axis = " << num << endl;
	imshow("symaxis image", img);

}
/*
void detectsymaxis(vector<Vec4i> in_lines)
{
	vector<Point> points;
	vector< vector<Point>> line(in_lines.size());
	for (int i = 0; i < in_lines.size(); i++)
		line[i].resize(in_lines[i].rows / 2);
	int n = 0;
	for (int i = 0; i < in_lines.size(); i++)
	{
		line[i][0].x = in_lines[n][0];
		line[i][0].y = in_lines[n][1];
		line[i][1].x = in_lines[n][2];
		line[i][1].y = in_lines[n][3];
		n++;
	}
	//ofstream fout("1.txt");
	vector<double> theta;
	for (int i = 0; i < line.size(); i++)
	{
		double dx = line[i][1].x - line[i][0].x;
		double dy = line[i][1].y - line[i][0].y;
		if (dx == 0)
		{
			theta.push_back(INFINITE_NUMBER);
			//fout << INFINITE_NUMBER << endl;
		}
		else if (dy == 0)
		{
			theta.push_back(0);
			//fout << 0 << endl;
		}
		else
		{
			theta.push_back(dy / dx);
			//fout << abs(dy / dx) << endl;
		}
	}
	//fout.close();
	for (int i = 0; i < theta.size() - 1; i++)
	{
		for (int j = 0; j < theta.size() - 1; j++)
		{
			if (theta[i] - theta[i + 1] < 0.02)
			{

			}
		}
	}
	int num = 0;
	for (int j = 0; j < theta.size() - 1; j++)
	{
		for (int i = 0; i < theta.size() - 1; i++)
		{
			if (abs(theta[j] - theta[i + 1]) < 0.000000001)
			{ 
				double b1 = line[j][1].y - theta[j] * line[j][1].x;
				double b2 = line[i + 1][1].y - theta[i + 1] * line[i + 1][1].x;
				double dis = abs(b1 - b2) / sqrt(theta[j] * theta[j] + 1);
				if (dis < 0.0001)
					continue;
				else
				{
					if (distance_equal)
						num++;
					else
						continue;
				}
				
			}
			else
				continue;
		}
	}
	
	cout << "Number of symmetry axis = " << num << endl;
	
}

static void on_trackbar_contours(int, void*)
{
    Mat cnt_img = Mat::zeros(w, w, CV_8UC3);
    int _levels = levels - 3;
    drawContours(cnt_img, contours, _levels <= 0 ? 3 : -1, Scalar(128, 255, 255),
                  3, LINE_AA, hierarchy, std::abs(_levels));

    imshow("Contours image", cnt_img);
}

void on_trackbar_binimg(int, void*)
{
	threshold(img, binimg, 100, 200, CV_THRESH_BINARY);
	imshow("二值化图像", binimg);
}
*/
int main()
{
	Mat img;
	Mat binimg;
	img = imread("111.jpg", 0);	//将读入的彩色图像直接以灰度图像读入
	namedWindow("原图", 1);
	imshow("原图", img);
	cout << "phase 1" << endl;
	binimg = img.clone();
	//进行二值化处理，选择119，200为阈值
	threshold(img, binimg, 119, 200, CV_THRESH_BINARY);
	namedWindow("二值化图像");
	imshow("二值化图像", binimg);
	cout << "phase 2" << endl;
	findContours(binimg, contours, hierarchy, RETR_TREE, CHAIN_APPROX_SIMPLE);
	binimg = cv::Scalar::all(0);
	drawContours(binimg, contours, -1, Scalar::all(255), 1, 8, hierarchy);
	imshow("Contours image", binimg);
	cout << "phase 3" << endl;
	//Canny边缘提取，
	Mat canny;
	Canny(binimg, canny, 350, 400);
	imshow("Canny边缘", canny);
	cout << "phase 4" << endl;
	//概率霍夫变换提取
	HoughLineFinder hlfinder;
	hlfinder.setLengthAndGap(10, 5);
	hlfinder.setminVote(5);
	vector<Vec4i> lines = hlfinder.findLines(canny);
	hlfinder.drawDetectedLines(binimg);
	imshow("HoughLine", binimg);
	cout << "phase 5" << endl;


	
	//createTrackbar("levels_bin", "二值化图像", 0, 255, on_trackbar_binimg);
	//on_trackbar_binimg(0, 0);
	
	

	//createTrackbar("levels_contours", "Contours image", &levels, 10, on_trackbar_contours);
	//on_trackbar_contours(0, 0);

	harris har;
	vector<Point> corners;
	goodFeaturesToTrack(binimg, corners, 20, 0.08, 10, noArray(), 12, false);
	har.drawOnImage(binimg, corners);
	imshow("Shi-Tomasi", binimg);
	cout << "phase 6" << endl;

	detectsymaxis(binimg, corners, binimg.rows, binimg.cols);
	waitKey();
	return 0;
}
