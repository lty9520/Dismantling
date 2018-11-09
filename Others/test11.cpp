#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>

using namespace std;

class mpoint{                       //class point(x, y, z ,G)
public:
	int x;
	int y;
	int z;
	double G;
	

};

void main()
{
	vector<mpoint> points;

	ifstream fin("222.txt");
	string  s;
	while (getline(fin, s))
	{
		vector<string> line;
		boost::split(line, s, boost::is_any_of("\t"), boost::token_compress_on);
		mpoint temp;
		temp.x = atoi(line[0].c_str());
		temp.y = atoi(line[1].c_str());
		temp.z = atoi(line[2].c_str());
		temp.G = atof(line[3].c_str());
		points.push_back(temp);
	}

	vector<double> G;
	vector<double> G_avr;
	int num = 0;

	mpoint p_pre;

	p_pre = points[0];
	G.push_back(points[0].G);

	for (int i = 1; i < points.size() - 2; i++)
	{
		if (abs(points[i].x - p_pre.x) >= 3 || 
			abs(points[i].y - p_pre.y) >= 2 ||
			abs(points[i].z - p_pre.z) >= 3)
		{
			double avr = 0;
			double temp = 0;
			for (int m = 0; m < num; m++)
			{
				temp = temp + G[m];
			}
			avr = temp / num;
			G_avr.push_back(avr);
			p_pre = points[i];
			G.clear();
			G.push_back(p_pre.G);
			num = 1;
		}
		else
		{
			G.push_back(points[i].G);
			num++;
		}
	}

	ofstream fout("G_avr.txt");
	fout << "id \t G_avr";
	fout << "\n";
	for (int i = 0; i < G_avr.size(); i++)
	{
		fout << i;
		fout << "\t";
		fout << G_avr[i];
		fout << "\n";
	}
	fout.close();

}