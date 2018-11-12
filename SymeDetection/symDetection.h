#ifndef SYMDETECTION_H
#define SYMDETECTION_H

#include <iostream>
#include <opencv2\highgui\highgui.hpp>
#include <opencv2\core\core.hpp>
#include <opencv2\imgproc\imgproc.hpp>
#include <math.h>
#include <fstream>

#define M_PI 3.14159265358979323846
#define INFINITE_NUMBER 10000

using namespace cv;
using namespace std;

class symetryDetection
{


public:


	//************************************
	// Method:    one_lines
	// FullName:  one_lines
	// Access:    public 
	// Returns:   bool
	// Qualifier:
	// Parameter: vector<Point> line1
	// Parameter: vector<Point> line2
	//************************************
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

	//************************************
	// Method:    distance_equal
	// FullName:  distance_equal
	// Access:    public 
	// Returns:   bool
	// Qualifier:
	// Parameter: vector<Point> line1
	// Parameter: vector<Point> line2
	//************************************
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

	//************************************
	// Method:    draw_symaxis
	// FullName:  draw_symaxis
	// Access:    public 
	// Returns:   void
	// Qualifier:
	// Parameter: Mat img
	// Parameter: Point l1_p1
	// Parameter: Point l1_p2
	// Parameter: Point l2_p1
	// Parameter: Point l2_p2
	//************************************
	void draw_symaxis(Mat img, Point l1_p1, Point l1_p2, Point l2_p1, Point l2_p2)
	{
		Point result_p1, result_p2;
		result_p1.x = (min_x(l1_p1, l1_p2) + min_x(l2_p1, l2_p2)) / 2;
		result_p1.y = (min_y(l1_p1, l1_p2) + min_y(l2_p1, l2_p2)) / 2;
		result_p2.x = (max_x(l1_p1, l1_p2) + max_x(l2_p1, l2_p2)) / 2;
		result_p2.y = (max_y(l1_p1, l1_p2) + max_y(l2_p1, l2_p2)) / 2;
		line(img, result_p1, result_p2, Scalar::all(125));
	}

	//Ѱ�����¶˵�����϶˵�
	vector<Point> get_boud_point(vector<Point> line1, vector<Point> line2)
	{
		vector<Point> result;
		return result;
	}

	//���Գ���
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
		cv::imshow("symaxis image", img);

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
	imshow("��ֵ��ͼ��", binimg);
	}
	*/


	//On Symmetry Detection
	//************************************
	// Method:    symmetrydetection
	// FullName:  symmetrydetection
	// Access:    public 
	// Returns:   void
	// Qualifier:
	// Parameter: vector<Point> con_r	�Ұ벿�߽��
	// Parameter: vector<Point> con_l	��벿�߽��
	// Parameter: int num_symax			���նԳ�������
	//************************************
	void symmetrydetection(vector<Point> con_r, vector<Point> con_l)
	{
		//ƽ������
		if (con_r.size() > con_l.size())
		{
			int dif_size = con_r.size() - con_l.size();
			int inter = floor(con_r.size() / dif_size);
			for (; dif_size > 0; dif_size--)
			{
				con_r.erase(con_r.begin() + dif_size * inter);
			}
		}
		else
		{
			int dif_size = con_l.size() - con_r.size();
			int inter = floor(con_l.size() / dif_size);
			for (; dif_size > 0; dif_size--)
			{
				con_l.erase(con_l.begin() + dif_size * inter);
			}
		}

		int num_symax = 0;
		vector<double> line_k;
		//�ж��Ƿ�Ϊ�Գ���
		for (int i = 0; i < con_l.size(); i++)
		{
			int n = 0;
			vector<double> line_para = get_line_para(con_l[i], con_r[i]);
			for (int j = 0, m = con_r.size() - 1; j < con_l.size(), m >= 0; j++, m--)
			{
				vector<double> dis_lr = dis_p2l(con_l[j], con_r[m], line_para[0], line_para[1]);
				if (dis_lr[0] == 0 || dis_lr[1] == 0)
					continue;
				if (flag_disequal(dis_lr[0], dis_lr[1]))
					n++;
				else
					continue;
			}
			if (n == con_r.size() - 5)
			{
				num_symax++;

				line_k.push_back(line_para[0]);
			}

		}
		//�����ܽ��ĶԳ���ϲ�
		for (int i = 0; i < line_k.size() - 1; i++)
		{
			double rad = get_lines_arctan(line_k[i], line_k[i + 1], 1);
			if (abs(rad) < 1)
			{
				num_symax--;
			}
		}
		cout << "Symmetry axis number :" << num_symax << endl;

	}

	/** �����ֵͼ�������
	* @param[in] img  ����Ĵ�����ͼ��
	* @param[out] center ��������
	* @retval 0  �����ɹ�
	* @retval -1 ����ʧ��
	* @note ����ͼ���Ƕ�ֵ��ͼ��
	* @note xc=M10/M00, yc=M01/M00, ���� Mx_order,y_order=SUMx,y(I(x,y)*x^x_order*y^y_order)
	*/
	void getGravityCenter(vector<vector<Point> > contours, Mat img, vector<Vec4i> hierarchy)
	{
		//���������� 	
		vector<Moments> mu(contours.size());
		for (int i = 0; i < contours.size(); i++)
		{
			mu[i] = moments(contours[i], false);
		}
		//�������������� 	
		vector<Point2f> mc(contours.size());
		for (int i = 0; i < contours.size(); i++)
		{
			mc[i] = Point2d(mu[i].m10 / mu[i].m00, mu[i].m01 / mu[i].m00);
		}
		for (int i = 0; i < contours.size(); i++)
		{
			Scalar color = Scalar(255, 0, 0);
			drawContours(img, contours, i, color, 2, 8, hierarchy, 0, Point());
			circle(img, mc[i], 5, Scalar(0, 0, 255), -1, 8, 0);
			rectangle(img, boundingRect(contours.at(i)), cvScalar(0, 255, 0));
			char tam[100];
			sprintf(tam, "(%0.0f,%0.0f)", mc[i].x, mc[i].y);
			putText(img, tam, Point(mc[i].x, mc[i].y), FONT_HERSHEY_SIMPLEX, 0.4, cvScalar(255, 0, 255), 1);
		}
		cv::imshow("center", img);
	}

	//����Сyֵ���
	//************************************
	// Method:    get_miny_point_id
	// FullName:  get_miny_point_id
	// Access:    public 
	// Returns:   int
	// Qualifier:
	// Parameter: mpoint * points	���е�����
	// Parameter: int size			�����С
	//************************************
	int get_miny_point_id(vector<Point> points, int size){ //get the point with min_y
		int i, min_id = 0;
		double miny = 10000;
		for (i = 0; i < size; i++){
			if (points[i].y < miny){
				miny = points[i].y;
				min_id = i;
			}
		}
		return min_id;
	}

	//������ֵ
	//************************************
	// Method:    get_cos
	// FullName:  get_cos
	// Access:    public 
	// Returns:   void
	// Qualifier:
	// Parameter: vector<Point>  points		���е�����
	// Parameter: vector<double> mcos		���е�����ֵ����
	// Parameter: int id			��Сyֵ����
	// Parameter: int size			�����С
	//************************************
	void get_cos(vector<Point> points, vector<double> mcos, int id, int size){  //get point's cos
		int i;
		double coss;
		for (i = 0; i < size; i++){
			if (i == id){
				mcos[i] = 2;
			}
			else{
				coss = (points[i].x - points[id].x) / sqrt((points[i].x - points[id].x) * (points[i].x - points[id].x) + (points[i].y - points[id].y) * (points[i].y - points[id].y));
				mcos[i] = coss;
			}
		}
	}

	//��λ�ǣ�rad��
	//************************************
	// Method:    calcAzimuthAngle
	// FullName:  calcAzimuthAngle
	// Access:    public 
	// Returns:   double
	// Qualifier:
	// Parameter: double x0		���ĵ�x
	// Parameter: double y0		���ĵ�y
	// Parameter: double x		�����x
	// Parameter: double y		�����y
	//************************************
	double calcAzimuthAngle(double x0, double y0, double x, double y)
	{
		double dx = x - x0;
		double dy = y - y0;
		//�����޽�
		double theta = atan((y - y0) / (x - x0));
		if (dx > 0 && dy > 0){			//��һ����
			return theta;
		}
		else if (dx < 0 && dy > 0){		//�ڶ�����
			return theta + M_PI;
		}
		else if (dx < 0 && dy < 0){		//��������
			return theta + M_PI;
		}
		else if (dx > 0 && dy < 0){		//��������
			return theta + (2 * M_PI);
		}
	}

	//��������ֵ�Ӵ�С����
	//************************************
	// Method:    sort_points_down
	// FullName:  sort_points_down
	// Access:    public 
	// Returns:   void
	// Qualifier:
	// Parameter: vector<Point>  points			���е�����
	// Parameter: vector<double> mcos			���е�����ֵ����
	// Parameter: int size				�����С
	//************************************
	void sort_points_down(vector<Point> points, vector<double> mcos, int size){   //sort the points
		int i, j;
		double temp_cos;
		Point temp_point;
		for (i = 0; i < size; i++){
			for (j = 0; j < size - i - 1; j++){      //bubble sorting
				if (mcos[j] < mcos[j + 1]){
					temp_cos = mcos[j];
					mcos[j] = mcos[j + 1];
					mcos[j + 1] = temp_cos;

					temp_point = points[j];
					points[j] = points[j + 1];
					points[j + 1] = temp_point;
				}
			}
		}
	}



	//************************************
	// Method:    get_pointlist_inner_contour2
	// FullName:  get_pointlist_inner_contour2
	// Access:    public 
	// Returns:   void
	// Qualifier:
	// Parameter: Mat src
	// Parameter: vector<Point> & contourlist
	//************************************
	void get_pointlist_inner_contour2(Mat src, vector<Point> &contourlist)
	{

		vector<vector<Point> > contours_out;
		vector<Vec4i> hierarchy_out;
		findContours(src, contours_out, hierarchy_out, RETR_EXTERNAL, CHAIN_APPROX_NONE);


		vector<vector<Point> > contours_all;
		vector<Vec4i> hierarchy_all;
		findContours(src, contours_all, hierarchy_all, RETR_TREE, CHAIN_APPROX_NONE);

		if (contours_all.size() == contours_out.size()) //û��������������ǰ����  
			cout << "without inner contours" << endl;

		for (int i = 0; i < contours_out.size(); i++)
		{
			int conloursize = contours_out[i].size();
			for (int j = 0; j < contours_all.size(); j++)
			{
				int tem_size = contours_all[j].size();
				if (conloursize == tem_size)
				{
					swap(contours_all[j], contours_all[contours_all.size() - 1]);
					contours_all.pop_back();
					break;
				}
			}
		}

		//contours_all��ֻʣ��������  
		//�����������  
		double maxarea = 0;
		int maxAreaIdx = 0;
		for (int index = contours_all.size() - 1; index >= 0; index--)
		{
			double tmparea = fabs(contourArea(contours_all[index]));
			if (tmparea > maxarea)
			{
				maxarea = tmparea;
				maxAreaIdx = index;//��¼���������������  
			}
		}

		contourlist = contours_all[maxAreaIdx];

	}

	//************************************
	// Method:    get_center
	// FullName:  get_center
	// Access:    public 
	// Returns:   Point
	// Qualifier:
	// Parameter: vector<Point> points
	//************************************
	Point get_center(vector<Point> points)
	{
		int ma_x, ma_y, mi_x, mi_y, size, x_cen, y_cen;
		size = points.size();
		vector<int>x(size);
		vector<int>y(size);
		for (int i = 0; i < size; i++)
		{
			x[i] = points[i].x;
			y[i] = points[i].y;
		}
		sort(x.begin(), x.end());
		sort(y.begin(), y.end());
		ma_x = x[size - 1];
		mi_x = x[0];
		ma_y = y[size - 1];
		mi_y = y[0];
		x_cen = (ma_x - mi_x) / 2;
		y_cen = (ma_y - mi_y) / 2;
		Point center;
		center.x = x_cen;
		center.y = y_cen;
		return center;
	}

	//coutours merge
	//************************************
	// Method:    coutour_merge
	// FullName:  coutour_merge
	// Access:    public 
	// Returns:   vector<vector<Point> >			���ձ߽������[[����], [�Ұ��]]
	// Qualifier:
	// Parameter: Mat img							����ͼ��
	// Parameter: vector<vector<Point> > contours	�߽��
	// Parameter: int r								ͼ��row
	//************************************
	vector<vector<Point> > coutour_merge(Mat img, vector<vector<Point> > contours, int r, vector<Vec4i> hierarchy)
	{

		//��ʾͼ��
		Mat img2 = img.clone();
		//��ϱ߽��洢hierarchy�㼶��ϵ
		vector<int>Index;
		for (int i = 0; i < contours.size() - 1;)
		{
			int next;
			if (i == 0 && hierarchy[i][0] == -1)
			{
				next = i + 1;
				i = next;
				continue;
			}
			if (hierarchy[i][0] != -1)
			{
				next = hierarchy[i][0];
				Index.push_back(i);
				i = next;
			}
		}

		//cout << contours[1][1] << endl;
		//�洢���б߽������
		int size = 0;
		//�߽������ͳ��
		for (int i = 0; i < Index.size(); i++)
		{
			size = size + contours[Index[i]].size();
		}
		//˳��߽��洢����
		vector<Point>contours_result(size);
		//˳��߽�����������
		int n = 0;
		//�������㼶�ı߽��ת��Ϊȥ���㼶��ϵ��˳��洢
		for (int i = 0; i < Index.size(); i++)
		{
			int num = contours[Index[i]].size();
			for (int j = 0; j < num; j++)
			{
				contours_result[n] = contours[Index[i]][j];
				n++;
			}
		}
		//�жϱ߽�������Ƿ���ż��
		if (contours_result.size() % 2 != 0)
		{
			contours_result.pop_back();
		}
		//�߽������
		Point center_point;
		//����߽���������������
		center_point = get_center(contours_result);
		/*
		//��ʼ��
		vector<Point> startpoint;

		for (int i = 0; i < contours_result.size(); i++)
		{
		if (contours_result[i].x == r / 2)
		{
		startpoint.push_back(contours_result[i]);
		}
		}
		*/

		/*
		//yֵ��С����
		int id = 0;
		//��yֵ��С����
		id = get_miny_point_id(contours_result, contours_result.size());
		//����cosֵ����
		vector<double> mcos(contours_result.size());
		//�������cosֵ
		get_cos(contours_result, mcos, id, contours_result.size());
		//��cosֵ���򣨴Ӵ�С��
		sort_points_down(contours_result, mcos, contours_result.size());
		*/

		//�洢���㷽λ�Ǻ͵������, [[Aziang][id]]
		vector<vector<double> > rad(contours_result.size(), vector<double>(2));
		//��λ�����������
		int num = 0;
		//�������������ĵ�ķ�λ��
		for (int i = 0; i < contours_result.size(); i++)
		{
			rad[i][0] = calcAzimuthAngle(center_point.x, center_point.y, contours_result[i].x, contours_result[i].y);
			rad[i][1] = num;
			num++;
		}
		//����λ�Ǵ�С�������򣨴�С����Ĭ�ϣ���
		sort(rad.begin(), rad.end());
		//�洢���ߵ���Ұ�ߵ������
		vector<Point> p_left, p_right_down, p_right_up;
		//���߽�㰴�շ�λ�ǻ���Ϊ���ߺ��Ұ�ߣ��ֱ������Ӧ����
		for (int i = 0; i < rad.size(); i++)
		{
			if (rad[i][0] <= M_PI * 3 / 2 && rad[i][0] >= M_PI / 2)
			{
				p_left.push_back(contours_result[rad[i][1]]);
			}
			else if (rad[i][0] < M_PI * 2 && rad[i][0] > M_PI * 3 / 2)
			{
				p_right_down.push_back(contours_result[rad[i][1]]);
			}
			else if (rad[i][0] < M_PI / 2 && rad[i][0] >= 0)
			{
				p_right_up.push_back(contours_result[rad[i][1]]);
			}
			else
				continue;
		}
		//���ߵ㷽λ�Ǵ洢������ʽͬ��
		vector<vector<double> > rad_left(p_left.size(), vector<double>(2));
		int num_left = 0;
		for (int i = 0; i < p_left.size(); i++)
		{
			rad_left[i][0] = calcAzimuthAngle(center_point.x, center_point.y, p_left[i].x, p_left[i].y);
			rad_left[i][1] = num_left;
			num_left++;
		}
		//�����ߵ㷽λ�ǽ�������
		sort(rad_left.begin(), rad_left.end());
		//�������߽߱��洢����
		vector<Point> contour_left;
		//��˳��洢���߽߱��
		for (int i = 0; i < rad_left.size(); i++)
		{
			contour_left.push_back(p_left[rad_left[i][1]]);
		}
		reverse(contour_left.begin(), contour_left.end());
		//�Ұ���°벿�㷽λ�Ǵ洢������ʽͬ��
		vector<vector<double> > rad_right_down(p_right_down.size(), vector<double>(2));
		int num_right_down = 0;
		for (int i = 0; i < p_right_down.size(); i++)
		{
			rad_right_down[i][0] = calcAzimuthAngle(center_point.x, center_point.y, p_right_down[i].x, p_right_down[i].y);
			rad_right_down[i][1] = num_right_down;
			num_right_down++;
		}
		//���Ұ���°벿�㷽λ�ǽ�������
		sort(rad_right_down.begin(), rad_right_down.end());
		//�����Ұ���°벿�߽��洢����
		vector<Point> contour_right_down;
		//��˳��洢�Ұ���°벿�߽��
		for (int i = 0; i < rad_right_down.size(); i++)
		{
			contour_right_down.push_back(p_right_down[rad_right_down[i][1]]);
		}
		//�Ұ���ϰ벿�㷽λ�Ǵ洢������ʽͬ��
		vector<vector<double> > rad_right_up(p_right_up.size(), vector<double>(2));
		int num_right_up = 0;
		for (int i = 0; i < p_right_up.size(); i++)
		{
			rad_right_up[i][0] = calcAzimuthAngle(center_point.x, center_point.y, p_right_up[i].x, p_right_up[i].y);
			rad_right_up[i][1] = num_right_up;
			num_right_up++;
		}
		//���Ұ���ϰ벿�㷽λ�ǽ�������
		sort(rad_right_up.begin(), rad_right_up.end());
		//�����Ұ���ϰ벿�߽��洢����
		vector<Point> contour_right_up;
		//��˳��洢�����ϰ벿�߽��
		for (int i = 0; i < rad_right_up.size(); i++)
		{
			contour_right_up.push_back(p_right_up[rad_right_up[i][1]]);
		}
		//reverse(contour_right_up.begin(), contour_right_up.end());
		//�����Ұ��ȫ���߽��洢����
		vector<Point> contour_right;
		for (int i = 0; i < contour_right_down.size(); i++)
		{
			contour_right.push_back(contour_right_down[i]);
		}
		for (int i = 0; i < contour_right_up.size(); i++)
		{
			contour_right.push_back(contour_right_up[i]);
		}
		reverse(contour_right.begin(), contour_right.end());
		//���Ʊ߽�㣬�Լ����
		for (int i = 0; i < contour_right.size(); i++)
		{
			circle(img, contour_right[i], 1, Scalar::all(125));
		}

		circle(img, contour_right[contour_right.size() - 1], 5, Scalar::all(255));

		imshow("right", img);
		//���Ʊ߽�㣬�Լ����
		for (int i = 0; i < contour_left.size(); i++)
		{
			circle(img2, contour_left[i], 1, Scalar::all(125));
		}

		circle(img2, contour_left[contour_left.size() - 1], 5, Scalar::all(255));

		imshow("left", img2);

		//���ձ߽������[[����], [�Ұ��]]
		vector<vector<Point> >result;
		result.push_back(contour_left);
		result.push_back(contour_right);

		return result;
	}

	//************************************
	// Method:    get_line_para
	// FullName:  symetryDetection::get_line_para
	// Access:    public 
	// Returns:   vector<double>	[k, b]
	// Qualifier:
	// Parameter: Point p1	ֱ�ߵ�һ��
	// Parameter: Point p2	ֱ�ߵڶ���
	//************************************
	vector<double> get_line_para(Point p1, Point p2)
	{
		double k, b;
		k = (p2.y - p1.y) / (p2.x - p1.x);
		b = p2.y - k * p2.x;
		vector<double> result;
		result.push_back(k);
		result.push_back(b);
		return result;
	}

	//************************************
	// Method:    dis_p2l
	// FullName:  symetryDetection::dis_p2l
	// Access:    public 
	// Returns:   vector<double>	[dis_l, dis_r]
	// Qualifier:
	// Parameter: Point p_l		ֱ����ߵĵ�
	// Parameter: Point p_r		ֱ���ұߵĵ�
	// Parameter: double k		ֱ��б��
	// Parameter: double b		ֱ�߽ؾ�
	//************************************
	vector<double> dis_p2l(Point p_l, Point p_r, double k, double b)
	{
		double dis_l = abs((p_l.y - k * p_l.x - b) / sqrt(1 + k * k));
		double dis_r = abs((p_r.y - k * p_r.x - b) / sqrt(1 + k * k));
		vector<double> result;
		result.push_back(dis_l);
		result.push_back(dis_r);
		return result;
	}

	//************************************
	// Method:    flag_disequal
	// FullName:  symetryDetection::flag_disequal
	// Access:    public 
	// Returns:   bool		���ξ����С����ֵ��true����֮��false
	// Qualifier:
	// Parameter: double dis1	��һ�ξ���
	// Parameter: double dis2	�ڶ��ξ���
	//************************************
	bool flag_disequal(double dis1, double dis2)
	{
		//cout << "dis = " << abs(dis1 - dis2) << endl;
		if (abs(dis1 - dis2) <= 15)
			return true;
		else
			return false;
	}
	//************************************
	// Method:    SkeletonTest
	// FullName:  symetryDetection::SkeletonTest
	// Access:    public 
	// Returns:   void
	// Qualifier:
	// Parameter: Mat & srcimage
	//************************************
	void SkeletonTest(Mat &srcimage)
	{
		int kernel[9];
		int nl = srcimage.rows;
		int nc = srcimage.cols;
		vector<Point> delete_list;
		int A, B;
		while (true)
		{
			for (int j = 1; j < nl - 1; j++)
			{
				uchar* data_pre = srcimage.ptr<uchar>(j - 1);
				uchar* data = srcimage.ptr<uchar>(j);
				uchar* data_next = srcimage.ptr<uchar>(j + 1);
				for (int i = 1; i < (nc - 1); i++)
				{
					if (data[i] == 255)
					{
						kernel[0] = 1;
						if (data_pre[i] == 255) kernel[1] = 1;
						else  kernel[1] = 0;
						if (data_pre[i + 1] == 255) kernel[2] = 1;
						else  kernel[2] = 0;
						if (data[i + 1] == 255) kernel[3] = 1;
						else  kernel[3] = 0;
						if (data_next[i + 1] == 255) kernel[4] = 1;
						else  kernel[4] = 0;
						if (data_next[i] == 255) kernel[5] = 1;
						else  kernel[5] = 0;
						if (data_next[i - 1] == 255) kernel[6] = 1;
						else  kernel[6] = 0;
						if (data[i - 1] == 255) kernel[7] = 1;
						else  kernel[7] = 0;
						if (data_pre[i - 1] == 255) kernel[8] = 1;
						else  kernel[8] = 0;

						B = 0;
						for (int k = 1; k < 9; k++)
						{
							B = B + kernel[k];
						}
						if ((B >= 2) && (B <= 6))
						{
							A = 0;
							if (!kernel[1] && kernel[2]) A++;
							if (!kernel[2] && kernel[3]) A++;
							if (!kernel[3] && kernel[4]) A++;
							if (!kernel[4] && kernel[5]) A++;
							if (!kernel[5] && kernel[6]) A++;
							if (!kernel[6] && kernel[7]) A++;
							if (!kernel[7] && kernel[8]) A++;
							if (!kernel[8] && kernel[1]) A++;
							//
							if (A == 1)
							{
								if ((kernel[1] * kernel[3] * kernel[5] == 0)
									&& (kernel[3] * kernel[5] * kernel[7] == 0))
								{
									delete_list.push_back(Point(i, j));
								}
							}
						}
					}
				}
			}
			int size = delete_list.size();
			if (size == 0)
			{
				break;
			}
			for (int n = 0; n < size; n++)
			{
				Point tem;
				tem = delete_list[n];
				uchar* data = srcimage.ptr<uchar>(tem.y);
				data[tem.x] = 0;
			}
			delete_list.clear();
			for (int j = 1; j < nl - 1; j++)
			{
				uchar* data_pre = srcimage.ptr<uchar>(j - 1);
				uchar* data = srcimage.ptr<uchar>(j);
				uchar* data_next = srcimage.ptr<uchar>(j + 1);
				for (int i = 1; i < (nc - 1); i++)
				{
					if (data[i] == 255)
					{
						kernel[0] = 1;
						if (data_pre[i] == 255) kernel[1] = 1;
						else  kernel[1] = 0;
						if (data_pre[i + 1] == 255) kernel[2] = 1;
						else  kernel[2] = 0;
						if (data[i + 1] == 255) kernel[3] = 1;
						else  kernel[3] = 0;
						if (data_next[i + 1] == 255) kernel[4] = 1;
						else  kernel[4] = 0;
						if (data_next[i] == 255) kernel[5] = 1;
						else  kernel[5] = 0;
						if (data_next[i - 1] == 255) kernel[6] = 1;
						else  kernel[6] = 0;
						if (data[i - 1] == 255) kernel[7] = 1;
						else  kernel[7] = 0;
						if (data_pre[i - 1] == 255) kernel[8] = 1;
						else  kernel[8] = 0;

						B = 0;
						for (int k = 1; k < 9; k++)
						{
							B = B + kernel[k];
						}
						if ((B >= 2) && (B <= 6))
						{
							A = 0;
							if (!kernel[1] && kernel[2]) A++;
							if (!kernel[2] && kernel[3]) A++;
							if (!kernel[3] && kernel[4]) A++;
							if (!kernel[4] && kernel[5]) A++;
							if (!kernel[5] && kernel[6]) A++;
							if (!kernel[6] && kernel[7]) A++;
							if (!kernel[7] && kernel[8]) A++;
							if (!kernel[8] && kernel[1]) A++;
							//
							if (A == 1)
							{
								if ((kernel[1] * kernel[3] * kernel[7] == 0)
									&& (kernel[1] * kernel[5] * kernel[7] == 0))
								{
									delete_list.push_back(Point(i, j));
								}
							}
						}
					}
				}
			}
			if (size == 0)
			{
				break;
			}
			for (int n = 0; n < size; n++)
			{
				Point tem;
				tem = delete_list[n];
				uchar* data = srcimage.ptr<uchar>(tem.y);
				data[tem.x] = 0;
			}
			delete_list.clear();
		}
		imshow("skeleton", srcimage);
	}

	//************************************
	// Method:    get_lines_arctan
	// FullName:  get_lines_arctan
	// Access:    public 
	// Returns:   double	������ֱ�߼�н�
	// Qualifier:
	// Parameter: float line_1_k	ֱ��1б��
	// Parameter: float line_2_k	ֱ��2б��
	// Parameter: int model			0Ϊ���ػ��ȣ�����Ϊ���ؽǶ�
	//************************************
	double get_lines_arctan(float line_1_k, float line_2_k, int model)
	{
		if (model == 0)
		{
			double tan_k = 0; //ֱ�߼н�����ֵ
			double lines_arctan;//ֱ��б�ʵķ�����ֵ
			tan_k = (line_2_k - line_1_k) / (1 + line_2_k*line_1_k); //��ֱ�߼нǵĹ�ʽ
			lines_arctan = atan(tan_k);
			return lines_arctan;
		}
		else
		{
			double tan_k = 0; //ֱ�߼н�����ֵ
			double lines_arctan;//ֱ��б�ʵķ�����ֵ
			tan_k = (line_2_k - line_1_k) / (1 + line_2_k*line_1_k); //��ֱ�߼нǵĹ�ʽ
			lines_arctan = atan(tan_k)* 180.0 / 3.1415926;

			return lines_arctan;
		}
	}

};

#endif	//SYMDETECTION_H