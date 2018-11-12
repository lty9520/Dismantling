#include "NeibAreaMarking.h"

nbaMarking::nbaMarking()
{

}

nbaMarking::~nbaMarking()
{

}

//************************************
// Method:    bwLabel
// FullName:  bwLabel
// Access:    public 
// Returns:   void
// Qualifier:
// Parameter: const Mat & imgBw
// Parameter: Mat & imgLabeled
//************************************
void nbaMarking::bwLabel(const Mat& imgBw, Mat & imgLabeled)
{
	// ��ͼ����Χ����һ��
	Mat imgClone = Mat(imgBw.rows + 1, imgBw.cols + 1, imgBw.type(), Scalar(0));
	imgBw.copyTo(imgClone(Rect(1, 1, imgBw.cols, imgBw.rows)));

	imgLabeled.create(imgClone.size(), imgClone.type());
	imgLabeled.setTo(Scalar::all(0));

	vector<vector<Point>> contours;
	vector<Vec4i> hierarchy;
	findContours(imgClone, contours, hierarchy, CV_RETR_CCOMP, CV_CHAIN_APPROX_NONE);

	vector<int> contoursLabel(contours.size(), 0);
	int numlab = 1;
	// �����Χ����
	for (vector<vector<Point>>::size_type i = 0; i < contours.size(); i++)
	{
		if (hierarchy[i][3] >= 0) // �и�����
		{
			continue;
		}
		for (vector<Point>::size_type k = 0; k != contours[i].size(); k++)
		{
			imgLabeled.at<uchar>(contours[i][k].y, contours[i][k].x) = numlab;
		}
		contoursLabel[i] = numlab++;
	}
	// ���������
	for (vector<vector<Point>>::size_type i = 0; i < contours.size(); i++)
	{
		if (hierarchy[i][3] < 0)
		{
			continue;
		}
		for (vector<Point>::size_type k = 0; k != contours[i].size(); k++)
		{
			imgLabeled.at<uchar>(contours[i][k].y, contours[i][k].x) = contoursLabel[hierarchy[i][3]];
		}
	}
	// ���������صı��
	for (int i = 0; i < imgLabeled.rows; i++)
	{
		for (int j = 0; j < imgLabeled.cols; j++)
		{
			if (imgClone.at<uchar>(i, j) != 0 && imgLabeled.at<uchar>(i, j) == 0)
			{
				imgLabeled.at<uchar>(i, j) = imgLabeled.at<uchar>(i, j - 1);
			}
		}
	}
	imgLabeled = imgLabeled(Rect(1, 1, imgBw.cols, imgBw.rows)).clone(); // ���߽�ü���1����
}

//************************************
// Method:    Seed_Filling
// FullName:  Seed_Filling
// Access:    public 
// Returns:   void
// Qualifier: //������䷨
// Parameter: const cv::Mat & binImg	��ֵͼ��
// Parameter: cv::Mat & lableImg		��Ǻ�ͼ��
// Parameter: int model_flag			8/4 ����
//************************************
void nbaMarking::Seed_Filling(const cv::Mat& binImg, cv::Mat& lableImg, int model_flag)   //������䷨
{

	if (model_flag != 8 && model_flag != 4)
	{
		cout << "Neighbour Area Model Input Error, Please Input Model As 8 or 4" << endl;
		return;
	}
	if (binImg.empty() ||
		binImg.type() != CV_8UC1)
	{
		return;
	}

	lableImg.release();
	binImg.convertTo(lableImg, CV_32SC1);

	//ofstream fout("5.txt");
	//for (int i = 0; i < lableImg.rows; i++)
	//{
	//	int* temp = lableImg.ptr<int>(i);
	//	for (int j = 0; j < lableImg.cols; j++)
	//	{
	//		int a = temp[j];
	//		fout << a;
	//		fout << "\t";
	//	}
	//	fout << "\n";
	//}
	//fout.close();

	int label = 1;

	int rows = binImg.rows - 1;
	int cols = binImg.cols - 1;
	if (model_flag == 8)
	{
		// 8�ڽӷ���
		for (int i = 1; i < rows - 1; i++)
		{
			int* data = lableImg.ptr<int>(i);
			for (int j = 1; j < cols - 1; j++)
			{
				if (data[j] == 255)
				{
					std::stack<std::pair<int, int>> neighborPixels;
					neighborPixels.push(std::pair<int, int>(i, j));     // ����λ��: <i,j>
					++label;  // û���ظ����ţ���ʼ�µı�ǩ
					while (!neighborPixels.empty())
					{
						std::pair<int, int> curPixel = neighborPixels.top(); //�������һ����һ�������غ���������һ�е��Ǹ��ŵı�Ÿ�����
						int curX = curPixel.first;
						int curY = curPixel.second;
						lableImg.at<int>(curX, curY) = label;

						neighborPixels.pop();

						//		��������������������8�����Ӧ��ϵ������������������
						//		|(X - 1, Y - 1)		(X - 1, Y)		(X - 1, Y + 1) |
						//		|  (X, Y - 1)		  (X, Y)		  (X, Y + 1)   |
						//		|(X + 1, Y - 1)		(X + 1, Y)		(X + 1, Y + 1) |
						//		����������������������������������������������������
						if (lableImg.at<int>(curX - 1, curY) == 255)
						{// �ϱ�
							neighborPixels.push(std::pair<int, int>(curX - 1, curY));
						}
						if (lableImg.at<int>(curX - 1, curY + 1) == 255)
						{// ���ϱ�
							neighborPixels.push(std::pair<int, int>(curX - 1, curY));
						}
						if (lableImg.at<int>(curX, curY + 1) == 255)
						{// �ұ�
							neighborPixels.push(std::pair<int, int>(curX, curY + 1));
						}
						if (lableImg.at<int>(curX + 1, curY + 1) == 255)
						{// ���±�
							neighborPixels.push(std::pair<int, int>(curX, curY + 1));
						}
						if (lableImg.at<int>(curX + 1, curY) == 255)
						{// �±�
							neighborPixels.push(std::pair<int, int>(curX + 1, curY));
						}
						if (lableImg.at<int>(curX + 1, curY - 1) == 255)
						{// ���±�
							neighborPixels.push(std::pair<int, int>(curX + 1, curY));
						}
						if (lableImg.at<int>(curX, curY - 1) == 255)
						{//���
							neighborPixels.push(std::pair<int, int>(curX, curY - 1));
						}
						if (lableImg.at<int>(curX - 1, curY - 1) == 255)
						{//���ϱ�
							neighborPixels.push(std::pair<int, int>(curX, curY - 1));
						}
					}
				}
			}
		}
	}
	else
	{
		// 4�ڽӷ���
		for (int i = 1; i < rows - 1; i++)
		{
			int* data = lableImg.ptr<int>(i);
			for (int j = 1; j < cols - 1; j++)
			{
				if (data[j] == 255)
				{
					std::stack<std::pair<int, int>> neighborPixels;
					neighborPixels.push(std::pair<int, int>(i, j));     // ����λ��: <i,j>
					++label;  // û���ظ����ţ���ʼ�µı�ǩ
					while (!neighborPixels.empty())
					{
						std::pair<int, int> curPixel = neighborPixels.top(); //�������һ����һ�������غ���������һ�е��Ǹ��ŵı�Ÿ�����
						int curX = curPixel.first;
						int curY = curPixel.second;
						lableImg.at<int>(curX, curY) = label;

						neighborPixels.pop();

						//		��������������������8�����Ӧ��ϵ������������������
						//		|(X - 1, Y - 1)		(X - 1, Y)		(X - 1, Y + 1) |
						//		|  (X, Y - 1)		  (X, Y)		  (X, Y + 1)   |
						//		|(X + 1, Y - 1)		(X + 1, Y)		(X + 1, Y + 1) |
						//		����������������������������������������������������
						if (lableImg.at<int>(curX - 1, curY) == 255)
						{// �ϱ�
							neighborPixels.push(std::pair<int, int>(curX - 1, curY));
						}
						if (lableImg.at<int>(curX, curY + 1) == 255)
						{// �ұ�
							neighborPixels.push(std::pair<int, int>(curX, curY + 1));
						}
						if (lableImg.at<int>(curX + 1, curY) == 255)
						{// �±�
							neighborPixels.push(std::pair<int, int>(curX + 1, curY));
						}
						if (lableImg.at<int>(curX, curY - 1) == 255)
						{//���
							neighborPixels.push(std::pair<int, int>(curX, curY - 1));
						}
					}
				}
			}
		}
	}



}


//************************************
// Method:    smallerone
// FullName:  smallerone
// Access:    public 
// Returns:   int
// Qualifier:
// Parameter: int num1
// Parameter: int num2
//************************************
int nbaMarking::smallerone(int num1, int num2)
{
	if (num1 < num2)
		return num1;
	else
		return num2;
}


//************************************
// Method:    Two_Pass
// FullName:  Two_Pass
// Access:    public 
// Returns:   void
// Qualifier: //����ɨ�跨
// Parameter: const cv::Mat & binImg
// Parameter: cv::Mat & lableImg
//************************************
void nbaMarking::Two_Pass(const cv::Mat& binImg, cv::Mat& lableImg)    //����ɨ�跨
{
	if (binImg.empty() ||
		binImg.type() != CV_8UC1)
	{
		return;
	}

	//IplImage binImgTemp = binImg;
	//src = cvCloneImage(&binImgTemp);
	//IplImage *pGrayImage_8U = cvCreateImage(cvGetSize(src), IPL_DEPTH_8U, 1);
	//IplImage *pGrayImage_32S = cvCreateImage(cvGetSize(src), IPL_DEPTH_32S, 1);
	//
	//cvConvertScale(pGrayImage_8U, pGrayImage_32S);  //8Uת32S
	// ��һ��ͨ·

	lableImg.release();
	binImg.convertTo(lableImg, CV_32SC1);

	//lableImg = cvarrToMat(pGrayImage_32S);

	int label = 1;
	std::vector<int> labelSet;
	labelSet.push_back(0);
	labelSet.push_back(1);

	//uchar *temp = lableImg.ptr<uchar>(26);//��ȡ�����е��׵�ַ
	//int a = temp[248];//��ȡ�����е����е�����ֵ����ֵ��a��
	//cout << a << endl;//�������ֵ
	//ofstream fout("1.txt");
	//
	//for (int nrow = 0; nrow < binImg.rows; nrow++)
	//{
	//	for (int ncol = 0; ncol < binImg.cols; ncol++)
	//	{
	//		uchar val = binImg.at<uchar>(nrow, ncol);
	//		if (val == 255)
	//		{
	//			val = 1;
	//		}
	//		else
	//		{
	//			val = 0;
	//		}
	//
	//		fout << int(val);		//������ֵת��Ϊ������д���ļ�
	//		//cout << int(val);
	//	}
	//	//cout << endl;			//��ӡһ�к���л���
	//	fout << endl;
	//}
	//
	//fout.close();

	//ofstream fout("3.txt");
	//for (int i = 0; i < lableImg.rows; i++)
	//{
	//	uchar* temp = lableImg.ptr<uchar>(i);
	//	for (int j = 0; j < lableImg.cols; j++)
	//	{
	//		int a = temp[j];
	//		fout << a;
	//		fout << "\t";
	//	}
	//	fout << "\n";
	//}
	//fout.close();


	int rows = binImg.rows - 1;
	int cols = binImg.cols - 1;
	for (int i = 1; i < rows; i++)
	{
		int* data_preRow = lableImg.ptr<int>(i - 1);
		int* data_curRow = lableImg.ptr<int>(i);
		for (int j = 1; j < cols; j++)
		{
			if (data_curRow[j] == 255)
			{
				std::vector<int> neighborLabels;
				neighborLabels.reserve(2);
				int leftPixel = data_curRow[j - 1];
				int upPixel = data_preRow[j];
				int leftupPixel = data_preRow[j - 1];
				int rightupPixel = data_preRow[j + 1];

				if (leftPixel > 1 && rightupPixel > 1)
				{
					if (leftPixel < rightupPixel)
					{
						neighborLabels.push_back(leftPixel);
						rightupPixel = leftPixel;
					}
					else
					{
						neighborLabels.push_back(rightupPixel);
						leftPixel = rightupPixel;
					}
				}
				if (leftupPixel > 1 && rightupPixel > 1)
				{
					if (leftupPixel < rightupPixel)
					{
						neighborLabels.push_back(leftupPixel);
						rightupPixel = leftupPixel;
					}
					else
					{
						neighborLabels.push_back(rightupPixel);
						leftupPixel = rightupPixel;
					}
				}

				if (neighborLabels.empty())
				{
					labelSet.push_back(++label);  // ����ͨ����ǩ+1
					data_curRow[j] = label;
					labelSet[label] = label;
				}
				else
				{
					if (leftPixel > 1)
					{
						neighborLabels.push_back(leftPixel);
					}
					else if (leftupPixel > 1)
					{
						neighborLabels.push_back(leftupPixel);
					}
					else if (upPixel > 1)
					{
						neighborLabels.push_back(upPixel);
					}
					else if (rightupPixel > 1)
					{
						neighborLabels.push_back(rightupPixel);
					}

					std::sort(neighborLabels.begin(), neighborLabels.end());
					int smallestLabel = neighborLabels[0];


					// ������С�ȼ۱�
					for (size_t k = 1; k < neighborLabels.size(); k++)
					{
						int tempLabel = neighborLabels[k];
						int& oldSmallestLabel = labelSet[tempLabel];
						if (oldSmallestLabel > smallestLabel)
						{
							labelSet[oldSmallestLabel] = smallestLabel;
							oldSmallestLabel = smallestLabel;
						}
						else if (oldSmallestLabel < smallestLabel)
						{
							labelSet[smallestLabel] = oldSmallestLabel;
						}
					}
				}
			}
		}
	}

	// ���µȼ۶��б�
	// ����С��Ÿ��ظ�����
	for (size_t i = 2; i < labelSet.size(); i++)
	{
		int curLabel = labelSet[i];
		int preLabel = labelSet[curLabel];
		while (preLabel != curLabel)
		{
			curLabel = preLabel;
			preLabel = labelSet[preLabel];
		}
		labelSet[i] = curLabel;
	};

	//uchar *temp1 = lableImg.ptr<uchar>(78);//��ȡ�����е��׵�ַ
	//int a1 = temp1[204];//��ȡ�����е����е�����ֵ����ֵ��a��
	//cout << a1 << endl;//�������ֵ

	//ofstream fout("3.txt");
	//for (int i = 0; i < lableImg.rows; i++)
	//{
	//	uchar* temp = lableImg.ptr<uchar>(i);
	//	for (int j = 0; j < lableImg.cols; j++)
	//	{
	//		int a = temp[j];
	//		fout << a;
	//		fout << "\t";
	//	}
	//	fout << "\n";
	//}
	//fout.close();


	for (int i = 0; i < rows; i++)
	{
		int* data = lableImg.ptr<int>(i);
		for (int j = 0; j < cols; j++)
		{
			int pixelLabel = data[j];
			pixelLabel = labelSet[pixelLabel];
		}
	}
}


/*

void Two_Pass(const cv::Mat& binImg, cv::Mat& lableImg)    //����ɨ�跨
{
if (binImg.empty() ||
binImg.type() != CV_8UC1)
{
return;
}

// ��һ��ͨ·

lableImg.release();
binImg.convertTo(lableImg, CV_32SC1);

int label = 1;
std::vector<int> labelSet;
labelSet.push_back(0);
labelSet.push_back(1);

int rows = binImg.rows - 1;
int cols = binImg.cols - 1;
for (int i = 1; i < rows; i++)
{
int* data_preRow = lableImg.ptr<int>(i - 1);
int* data_curRow = lableImg.ptr<int>(i);
for (int j = 1; j < cols; j++)
{
if (data_curRow[j] == 1)
{
std::vector<int> neighborLabels;
neighborLabels.reserve(2);
int leftPixel = data_curRow[j - 1];
int upPixel = data_preRow[j];
if (leftPixel > 1)
{
neighborLabels.push_back(leftPixel);
}
if (upPixel > 1)
{
neighborLabels.push_back(upPixel);
}

if (neighborLabels.empty())
{
labelSet.push_back(++label);  // ����ͨ����ǩ+1
data_curRow[j] = label;
labelSet[label] = label;
}
else
{
std::sort(neighborLabels.begin(), neighborLabels.end());
int smallestLabel = neighborLabels[0];
data_curRow[j] = smallestLabel;

// ������С�ȼ۱�
for (size_t k = 1; k < neighborLabels.size(); k++)
{
int tempLabel = neighborLabels[k];
int& oldSmallestLabel = labelSet[tempLabel];
if (oldSmallestLabel > smallestLabel)
{
labelSet[oldSmallestLabel] = smallestLabel;
oldSmallestLabel = smallestLabel;
}
else if (oldSmallestLabel < smallestLabel)
{
labelSet[smallestLabel] = oldSmallestLabel;
}
}
}
}
}
}

// ���µȼ۶��б�
// ����С��Ÿ��ظ�����
for (size_t i = 2; i < labelSet.size(); i++)
{
int curLabel = labelSet[i];
int preLabel = labelSet[curLabel];
while (preLabel != curLabel)
{
curLabel = preLabel;
preLabel = labelSet[preLabel];
}
labelSet[i] = curLabel;
};

for (int i = 0; i < rows; i++)
{
int* data = lableImg.ptr<int>(i);
for (int j = 0; j < cols; j++)
{
int& pixelLabel = data[j];
pixelLabel = labelSet[pixelLabel];
}
}
}

*/

//��ɫ��ʾ
//************************************
// Method:    GetRandomColor
// FullName:  GetRandomColor
// Access:    public 
// Returns:   cv::Scalar
// Qualifier:
//************************************
cv::Scalar GetRandomColor()
{
	uchar r = 255 * (rand() / (1.0 + RAND_MAX));
	uchar g = 255 * (rand() / (1.0 + RAND_MAX));
	uchar b = 255 * (rand() / (1.0 + RAND_MAX));
	return cv::Scalar(b, g, r);
}


//************************************
// Method:    LabelColor
// FullName:  LabelColor
// Access:    public 
// Returns:   void
// Qualifier:
// Parameter: const cv::Mat & labelImg
// Parameter: cv::Mat & colorLabelImg
//************************************
void nbaMarking::LabelColor(const cv::Mat& labelImg, cv::Mat& colorLabelImg)
{
	if (labelImg.empty() ||
		labelImg.type() != CV_32SC1
		//labelImg.type() != CV_8UC1
		)
	{
		return;
	}

	std::map<int, cv::Scalar> colors;

	int rows = labelImg.rows;
	int cols = labelImg.cols;

	colorLabelImg.release();
	colorLabelImg.create(rows, cols, CV_8UC3);
	colorLabelImg = cv::Scalar::all(0);

	for (int i = 0; i < rows; i++)
	{
		const int* data_src = (int*)labelImg.ptr<int>(i);
		uchar* data_dst = colorLabelImg.ptr<uchar>(i);
		for (int j = 0; j < cols; j++)
		{
			int pixelValue = data_src[j];
			if (pixelValue > 1)
			{
				if (colors.count(pixelValue) <= 0)
				{
					colors[pixelValue] = GetRandomColor();
				}

				cv::Scalar color = colors[pixelValue];
				*data_dst++ = color[0];
				*data_dst++ = color[1];
				*data_dst++ = color[2];
			}
			else
			{
				data_dst++;
				data_dst++;
				data_dst++;
			}
		}
	}
}

