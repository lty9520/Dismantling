#include "..\Tools\benchDetect.h"


benchDetect::benchDetect()
{

}

benchDetect::~benchDetect()
{

}



//分支点标记
//************************************
// Method:    bench_delete
// FullName:  bench_delete
// Access:    public 
// Returns:   void
// Qualifier:
// Parameter: Mat srcimg
//************************************
void benchDetect::bench_detect(Mat srcimg)
{
	//骨架端点向量
	vector<Point> end_point;
	//骨架分支点向量
	vector<Point> bench_point;
	//既不是分支点也不是端点向量
	vector<Point> other_point;
	//正常骨架连续点向量
	vector<Point> skele_point;
	int neighbourhood[9];
	int row = srcimg.rows;
	int col = srcimg.cols;
	int num = 0;
	for (int j = 1; j < row - 1; j++)
	{
		uchar* data_last = srcimg.ptr<uchar>(j - 1);
		uchar* data = srcimg.ptr<uchar>(j);
		uchar* data_next = srcimg.ptr<uchar>(j + 1);
		for (int i = 1; i < col - 1; i++)
		{
			if (data[i] == 255)
			{
				//	———————————————————————对应关系———————————————————————
				//	|p点8邻域图	|	neighbourhood[]		|								data[]							|
				//	|9	2	3	|		8	1	2		|	data_last[i - 1]		data_last[i]		data_last[i + 1]|
				//	|8	1	4	|		7	0	3		|	  data[i - 1]				data[i]				data[i + 1]	|
				//	|7	6	5	|		6	5	4		|	data_next[i - 1]		data_next[i]		data_next[i + 1]|
				//  ——————————————————————————————————————————————————
				neighbourhood[0] = 1;

				if (data_last[i] == 255)
					neighbourhood[1] = 1;
				else
					neighbourhood[1] = 0;

				if (data_last[i + 1] == 255)
					neighbourhood[2] = 1;
				else
					neighbourhood[2] = 0;

				if (data[i + 1] == 255)
					neighbourhood[3] = 1;
				else
					neighbourhood[3] = 0;

				if (data_next[i + 1] == 255)
					neighbourhood[4] = 1;
				else
					neighbourhood[4] = 0;

				if (data_next[i] == 255)
					neighbourhood[5] = 1;
				else
					neighbourhood[5] = 0;

				if (data_next[i - 1] == 255)
					neighbourhood[6] = 1;
				else
					neighbourhood[6] = 0;

				if (data[i - 1] == 255)
					neighbourhood[7] = 1;
				else
					neighbourhood[7] = 0;

				if (data_last[i - 1] == 255)
					neighbourhood[8] = 1;
				else
					neighbourhood[8] = 0;

				int neighbour_sum = 0;
				for (int k = 1; k < 9; k++)
				{
					neighbour_sum += neighbourhood[k];
				}

				if (neighbourhood[0] == 1 && neighbour_sum < 2)
					end_point.push_back(Point(i, j));
				else if (neighbourhood[0] == 1 && neighbour_sum > 2)
				{
					bench_point.push_back(Point(i, j));
					data[i] = 0;
					num++;
				}
				else if (neighbourhood[0] == 1 && neighbour_sum == 2 &&
					(neighbourhood[1] + neighbourhood[2] == 2 ||
					neighbourhood[2] + neighbourhood[3] == 2 ||
					neighbourhood[3] + neighbourhood[4] == 2 ||
					neighbourhood[4] + neighbourhood[5] == 2 ||
					neighbourhood[5] + neighbourhood[6] == 2 ||
					neighbourhood[6] + neighbourhood[7] == 2 ||
					neighbourhood[7] + neighbourhood[8] == 2
					))
					other_point.push_back(Point(i, j));
				else
					skele_point.push_back(Point(i, j));
			}
		}
	}
	cout << "分支点个数:" << num << endl;

}

