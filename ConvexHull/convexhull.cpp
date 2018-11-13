#include "convexhull.h"

getConvexhull::getConvexhull()
{

}

getConvexhull::~getConvexhull()
{

}

//求最小y值点号
//************************************
// Method:    get_miny_point_id
// FullName:  get_miny_point_id
// Access:    public 
// Returns:   int
// Qualifier:
// Parameter: mpoint * points	所有点数组
// Parameter: int size			数组大小
//************************************
int getConvexhull::get_miny_point_id(mpoint *points, int size){ //get the point with min_y
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

//求余弦值
//************************************
// Method:    get_cos
// FullName:  get_cos
// Access:    public 
// Returns:   void
// Qualifier:
// Parameter: mpoint * points	所有点数组
// Parameter: double * mcos		所有点余弦值数组
// Parameter: int id			最小y值点点号
// Parameter: int size			数组大小
//************************************
void getConvexhull::get_cos(mpoint *points, double *mcos, int id, int size){  //get point's cos
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

//求方位角（rad）
//************************************
// Method:    calcAzimuthAngle
// FullName:  calcAzimuthAngle
// Access:    public 
// Returns:   double
// Qualifier:
// Parameter: double x0		中心点x
// Parameter: double y0		中心点y
// Parameter: double x		所求点x
// Parameter: double y		所求点y
//************************************
double getConvexhull::calcAzimuthAngle(double x0, double y0, double x, double y)
{
	double dx = x - x0;
	double dy = y - y0;
	//求象限角
	double theta = atan((y - y0) / (x - x0));
	if (dx > 0 && dy > 0){			//第一象限
		return theta;
	}
	else if (dx < 0 && dy > 0){		//第二象限
		return theta + M_PI;
	}
	else if (dx < 0 && dy < 0){		//第三象限
		return theta + M_PI;
	}
	else if (dx > 0 && dy < 0){		//第四象限
		return theta + (2 * M_PI);
	}
}

//按照余弦值从大到小排序
//************************************
// Method:    sort_points_down
// FullName:  sort_points_down
// Access:    public 
// Returns:   void
// Qualifier:
// Parameter: mpoint * points		所有点数组
// Parameter: double * mcos			所有点余弦值数组
// Parameter: int size				数组大小
//************************************
void getConvexhull::sort_points_down(mpoint *points, double *mcos, int size){   //sort the points
	int i, j;
	double temp_cos;
	mpoint temp_point;
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

//返回逆时针方向的点
//************************************
// Method:    ccw
// FullName:  ccw
// Access:    public 
// Returns:   int
// Qualifier:
// Parameter: mpoint a		起始点
// Parameter: mpoint b		中间点
// Parameter: mpoint c		末尾点
//************************************
int getConvexhull::ccw(mpoint a, mpoint b, mpoint c){          //judge if it is couter-colockwise
	double area2 = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
	if (area2 < 0){
		return -1;          // clockwise
	}
	else{
		if (area2 > 0) return 1;    // counter-clockwise
		else return 0;              // collinear
	}

}

//按照余弦值从小到大排序
//************************************
// Method:    sort_points_up
// FullName:  sort_points_up
// Access:    public 
// Returns:   void
// Qualifier:
// Parameter: mpoint * points		所有点数组
// Parameter: double * mcos			所有点余弦值数组
// Parameter: int size				数组大小
//************************************
void getConvexhull::sort_points_up(mpoint *points, double *mcos, int size){   //sort the points
	int i, j;
	double temp_cos;
	mpoint temp_point;
	for (i = 0; i < size; i++){
		for (j = 0; j < size - i - 1; j++){      //bubble sorting
			if (mcos[j] > mcos[j + 1]){
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

//返回顺时针方向的点
//************************************
// Method:    cw
// FullName:  cw
// Access:    public 
// Returns:   int
// Qualifier:
// Parameter: mpoint a		起始点
// Parameter: mpoint b		中间点
// Parameter: mpoint c		末尾点
//************************************
int getConvexhull::cw(mpoint a, mpoint b, mpoint c){          //judge if it is couter-colockwise
	double area2 = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
	if (area2 > 0){
		return -1;          // counter-clockwise
	}
	else{
		if (area2 < 0) return 1;    // clockwise
		else return 0;              // collinear
	}

}

//得到ccw的结果点
//************************************
// Method:    get_outpoint_down
// FullName:  get_outpoint_down
// Access:    public 
// Returns:   void
// Qualifier:
// Parameter: mpoint * points			所有点数组
// Parameter: int size					数组大小
// Parameter: pcl::PointCloud<pcl::PointXYZ>::Ptr & cloud_tubao		结果点云
//************************************
void getConvexhull::get_outpoint_down(mpoint *points, int size, pcl::PointCloud<pcl::PointXYZ>::Ptr& cloud_tubao){    //get points in stack
	int i, k;
	vector <mpoint>outpoint;
	outpoint.push_back(points[0]);
	outpoint.push_back(points[1]);
	i = 2;
	while (true){
		if (i == size){
			break;
		}
		if (ccw(outpoint[outpoint.size() - 2], outpoint[outpoint.size() - 1], points[i]) > 0){
			outpoint.push_back(points[i]);
			i = i + 1;
		}
		else if (ccw(outpoint[outpoint.size() - 2], outpoint[outpoint.size() - 1], points[i]) == 0)
		{
			i++;
			continue;
		}
		else{
			outpoint.pop_back();
		}
	}
	cloud_tubao->width = outpoint.size();
	cloud_tubao->height = 1;
	cloud_tubao->points.resize(cloud_tubao->width * cloud_tubao->height);
	cloud_tubao->is_dense = true;

	cout << "The outpoints have: " << outpoint.size() << endl;
	for (k = 0; k < outpoint.size(); k++){
		cout << outpoint[k].x << " " << outpoint[k].y << endl;
		cloud_tubao->points[k].x = outpoint[k].x;
		cloud_tubao->points[k].y = outpoint[k].y;
	}


}

//得到cw的结果点
//************************************
// Method:    get_outpoint_up
// FullName:  get_outpoint_up
// Access:    public 
// Returns:   void
// Qualifier:	
// Parameter: mpoint * points		所有点数组
// Parameter: int size				数组大小
// Parameter: pcl::PointCloud<pcl::PointXYZ>::Ptr & cloud_tubao		结果点云
//************************************
void getConvexhull::get_outpoint_up(mpoint *points, int size, pcl::PointCloud<pcl::PointXYZ>::Ptr& cloud_tubao){    //get points in stack
	int i, k;
	vector <mpoint>outpoint;
	outpoint.push_back(points[0]);
	outpoint.push_back(points[1]);
	i = 2;
	while (true){
		if (i == size){
			break;
		}
		if (cw(outpoint[outpoint.size() - 2], outpoint[outpoint.size() - 1], points[i]) < 0){
			outpoint.push_back(points[i]);
			i = i + 1;
		}
		else if (cw(outpoint[outpoint.size() - 2], outpoint[outpoint.size() - 1], points[i]) == 0)
		{
			i++;
			continue;
		}
		else{
			outpoint.pop_back();
		}
		if (outpoint.size() <= 1)
		{
			break;
		}
	}
	cloud_tubao->width = outpoint.size();
	cloud_tubao->height = 1;
	cloud_tubao->points.resize(cloud_tubao->width * cloud_tubao->height);
	cloud_tubao->is_dense = true;

	cout << "The outpoints have: " << outpoint.size() << endl;
	for (k = 0; k < outpoint.size(); k++){
		//cout << outpoint[k].x << " " << outpoint[k].y << endl;
		cloud_tubao->points[k].x = outpoint[k].x;
		cloud_tubao->points[k].y = outpoint[k].y;
	}


}


//求数组最小值
//************************************
// Method:    min_array
// FullName:  min_array
// Access:    public 
// Returns:   double
// Qualifier: 
// Parameter: double * array	所要求的数组
// Parameter: int n				元素个数
//************************************
double
getConvexhull::min_array(double *array, int n)
{
	double min = array[0];
	for (int i = 0; i < n; i++)
	{
		if (min > array[i])
		{
			min = array[i];
		}
	}
	return min;
}

//求数组最大值
//************************************
// Method:    max_array
// FullName:  max_array
// Access:    public 
// Returns:   double
// Qualifier:
// Parameter: double * array	所要求的数组
// Parameter: int n				元素个数
//************************************
double
getConvexhull::max_array(double *array, int n)
{
	double max = array[0];
	for (int i = 0; i < n; i++)
	{
		if (max < array[i])
		{
			max = array[i];
		}
	}
	return max;
}

//求两点距离
//************************************
// Method:    getDistance
// FullName:  getDistance
// Access:    public 
// Returns:   double
// Qualifier:
// Parameter: double x1		第一点x
// Parameter: double y1		第一点y
// Parameter: double x2		第二点x
// Parameter: double y2		第二点y
//************************************
double getConvexhull::getDistance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}

//显示点号
//************************************
// Method:    showID
// FullName:  showID
// Access:    public 
// Returns:   void
// Qualifier:
// Parameter: pcl::PointCloud<pcl::PointXYZ>::Ptr & cloud			显示点云
// Parameter: boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer	显示对象
// Parameter: double scale		字号
// Parameter: double r			颜色r
// Parameter: double g			颜色g
// Parameter: double b			颜色b
//************************************
void getConvexhull::showID(pcl::PointCloud<pcl::PointXYZ>::Ptr& cloud,
	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer,
	double scale,
	double r,
	double g,
	double b)
{
	for (int i = 0; i < cloud->points.size() - 1; i++)
		viewer->addText3D(to_string(i), cloud->points[i], scale, r, g, b);
}

//画顺序点的连线
//************************************
// Method:    drawLine
// FullName:  drawLine
// Access:    public 
// Returns:   void
// Qualifier:
// Parameter: pcl::PointCloud<pcl::PointXYZ>::Ptr & cloud		显示点云
// Parameter: boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer		显示对象
// Parameter: double r			颜色r
// Parameter: double g			颜色g
// Parameter: double b			颜色b
//************************************
void getConvexhull::drawLine(pcl::PointCloud<pcl::PointXYZ>::Ptr& cloud,
	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer,
	double r,
	double g,
	double b)
{
	for (int i = 0; i < cloud->points.size() - 1; i++)
		viewer->addLine<pcl::PointXYZ>(cloud->points[i], cloud->points[i + 1], r, g, b, to_string(i));
	viewer->addLine<pcl::PointXYZ>(cloud->points[cloud->points.size() - 1], cloud->points[0], r, g, b, to_string(cloud->points.size() + 1));
}

//************************************
// Method:    calcConvexHull
// FullName:  calcConvexHull
// Access:    public 
// Returns:   void
// Qualifier:
// Parameter: double * x	投影方向x轴
// Parameter: double * y	投影方向y轴
// Parameter: int size		点云大小
// Parameter: pcl::PointCloud<pcl::PointXYZ>::Ptr & cloud_final	最终结果点云
//************************************
void getConvexhull::calcConvexHull(double *x,
	double *y,
	int size,
	pcl::PointCloud<pcl::PointXYZ>::Ptr& cloud_final)
{
	cloud_final->clear();
	//ccw计算结果点云
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_tubao_ccw(new pcl::PointCloud<pcl::PointXYZ>);
	cloud_tubao_ccw->clear();
	//cw计算结果点云
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_tubao_cw(new pcl::PointCloud<pcl::PointXYZ>);
	cloud_tubao_cw->clear();
	//整体结果点云
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_tubao_result(new pcl::PointCloud<pcl::PointXYZ>);
	cloud_tubao_result->clear();
	//点云平面中心点
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_center(new pcl::PointCloud<pcl::PointXYZ>);
	cloud_center->clear();

	cloud_center->width = 1;
	cloud_center->height = 1;
	cloud_center->points.resize(cloud_center->width * cloud_center->height);
	cloud_center->is_dense = true;


	cloud_tubao_ccw->width = size;
	cloud_tubao_ccw->height = 1;
	cloud_tubao_ccw->points.resize(cloud_tubao_ccw->width * cloud_tubao_ccw->height);
	cloud_tubao_ccw->is_dense = true;

	cloud_tubao_cw->width = size;
	cloud_tubao_cw->height = 1;
	cloud_tubao_cw->points.resize(cloud_tubao_cw->width * cloud_tubao_cw->height);
	cloud_tubao_cw->is_dense = true;

	cloud_tubao_result->width = size * 2;
	cloud_tubao_result->height = 1;
	cloud_tubao_result->points.resize(cloud_tubao_result->width * cloud_tubao_result->height);
	cloud_tubao_result->is_dense = true;
	//求x的最大值
	double max_x = max_array(x, size);
	//求x的最小值
	double min_x = min_array(x, size);
	//求y的最大值
	double max_y = max_array(y, size);
	//求y的最小值
	double min_y = min_array(y, size);
	//求x中间值
	double mid_x = (abs(max_x) - abs(min_x)) / 2;
	//求y中间值
	double mid_y = (abs(max_y) - abs(min_y)) / 2;

	cloud_center->points[0].x = mid_x;
	cloud_center->points[0].y = mid_y;

	for (int i = 0; i < size; i++)
	{
		cloud_tubao_ccw->points[i].x = x[i];
		cloud_tubao_ccw->points[i].y = y[i];
		cloud_tubao_ccw->points[i].z = 0;
	}

	for (int i = 0; i < size; i++)
	{
		cloud_tubao_cw->points[i].x = x[i];
		cloud_tubao_cw->points[i].y = y[i];
		cloud_tubao_cw->points[i].z = 0;
	}
	//计算ch的数组
	mpoint *points;
	points = new mpoint[size];
	for (int i = 0; i < size; i++)
	{
		points[i].x = x[i];
		points[i].y = y[i];
	}
	//**************ccw***************
	int miny_point_id_ccw;
	double *mcos_ccw;
	mcos_ccw = new double[size];
	//求最低点（ccw）
	miny_point_id_ccw = get_miny_point_id(points, size);
	//求cos值（ccw）
	get_cos(points, mcos_ccw, miny_point_id_ccw, size);
	//将点按cos排序（ccw）
	sort_points_down(points, mcos_ccw, size);
	//计算ch（ccw）
	get_outpoint_down(points, size, cloud_tubao_ccw);
	//*************cw*************
	int miny_point_id_cw;
	double *mcos_cw;
	mcos_cw = new double[size];
	//求最低点（cw）
	miny_point_id_cw = get_miny_point_id(points, size);
	//求cos值（cw）
	get_cos(points, mcos_cw, miny_point_id_cw, size);
	//将点按cos排序（cw）
	sort_points_up(points, mcos_cw, size);
	//计算ch（cw）
	get_outpoint_up(points, size, cloud_tubao_cw);
	//resize整体结果点云
	cloud_tubao_result->width = (cloud_tubao_cw->points.size() + cloud_tubao_ccw->points.size());
	cloud_tubao_result->height = 1;
	cloud_tubao_result->points.resize(cloud_tubao_result->width * cloud_tubao_result->height);
	cloud_tubao_result->is_dense = true;

	cloud_final->width = cloud_tubao_result->points.size();
	cloud_final->height = 1;
	cloud_final->points.resize(cloud_final->width * cloud_final->height);
	cloud_final->is_dense = true;
	//ccw结果转存
	for (int i = 0; i < cloud_tubao_cw->points.size(); i++)
	{
		cloud_tubao_result->points[i].x = cloud_tubao_cw->points[i].x;
		cloud_tubao_result->points[i].y = cloud_tubao_cw->points[i].y;
		cloud_tubao_result->points[i].z = cloud_tubao_cw->points[i].z;
	}
	//cw结果转存

	int j = 0;
	for (int i = cloud_tubao_cw->points.size(); i < (cloud_tubao_cw->points.size() + cloud_tubao_ccw->points.size()); i++)
	{
		cloud_tubao_result->points[i].x = cloud_tubao_ccw->points[j].x;
		cloud_tubao_result->points[i].y = cloud_tubao_ccw->points[j].y;
		cloud_tubao_result->points[i].z = cloud_tubao_ccw->points[j].z;
		j++;
	}
	//方位角存贮数组[[aa], [id]]
	vector<vector<double> > rad(cloud_tubao_result->points.size(), vector<double>(2));
	int num = 0;
	for (int i = 0; i < cloud_tubao_result->points.size(); i++){
		rad[i][0] = calcAzimuthAngle(mid_x, mid_y, cloud_tubao_result->points[i].x, cloud_tubao_result->points[i].y);
		rad[i][1] = num;
		num++;
	}
	//方位角排序
	sort(rad.begin(), rad.end());
	//最终ch结果转存


	for (int i = 0; i < cloud_final->points.size(); i++)
	{
		cloud_final->points[i].x = cloud_tubao_result->points[rad[i][1]].x;
		cloud_final->points[i].y = cloud_tubao_result->points[rad[i][1]].y;
		cloud_final->points[i].z = cloud_tubao_result->points[rad[i][1]].z;
	}

}

//************************************
// Method:    ch_operation
// FullName:  ch_operation
// Access:    public 
// Returns:   void
// Qualifier:
// Parameter: string file	文件名称(路径)
//************************************
void getConvexhull::ch_operation(string file)
{
	DWORD start_time = GetTickCount(); //开始计时

	//初始化PointCloud对象
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::io::loadPCDFile(file, *cloud);
	int size = cloud->points.size();
	//std::cout << "PointCloud before filtering has: " << cloud->points.size() << " data points." << std::endl;
	std::cerr << "Point cloud data: " << cloud->points.size() << " points" << std::endl;

	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_proj(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_tubao_ccw(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_tubao_cw(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_tubao_result(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_tubao_final(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_center(new pcl::PointCloud<pcl::PointXYZ>);

	cloud_center->width = 1;
	cloud_center->height = 1;
	cloud_center->points.resize(cloud_center->width * cloud_center->height);
	cloud_center->is_dense = true;


	cloud_tubao_ccw->width = size;
	cloud_tubao_ccw->height = 1;
	cloud_tubao_ccw->points.resize(cloud_tubao_ccw->width * cloud_tubao_ccw->height);
	cloud_tubao_ccw->is_dense = true;

	cloud_tubao_cw->width = size;
	cloud_tubao_cw->height = 1;
	cloud_tubao_cw->points.resize(cloud_tubao_cw->width * cloud_tubao_cw->height);
	cloud_tubao_cw->is_dense = true;

	cloud_proj->width = size;
	cloud_proj->height = 1;
	cloud_proj->points.resize(cloud_proj->width * cloud_proj->height);
	cloud_proj->is_dense = true;

	cloud_tubao_result->width = size * 2;
	cloud_tubao_result->height = 1;
	cloud_tubao_result->points.resize(cloud_tubao_result->width * cloud_tubao_result->height);
	cloud_tubao_result->is_dense = true;

	cout << " 0" << endl;
	double* z = new double[size];
	double* x = new double[size];
	double* y = new double[size];
	for (int i = 0; i < size; i++){

		x[i] = cloud->points[i].x;
		y[i] = cloud->points[i].y;
		z[i] = cloud->points[i].z;
		//fout << x[i];
		//fout << "\n";
		//fout << flush;
	}


	//求x的最大值
	double max_x = max_array(x, size);
	//求x的最小值
	double min_x = min_array(x, size);
	//求y的最大值
	double max_y = max_array(y, size);
	//求y的最小值
	double min_y = min_array(y, size);
	//求z的最大值
	double max_z = max_array(z, size);
	//求z的最小值
	double min_z = min_array(z, size);
	//求x中间值
	double mid_x = (max_x + min_x) / 2;
	//求y中间值			    
	double mid_y = (max_y + min_y) / 2;
	//求z中间值			    
	double mid_z = (max_z + min_z) / 2;

	for (int i = 0; i < cloud->points.size(); i++)
	{
		cloud_proj->points[i].x = cloud->points[i].y;
		cloud_proj->points[i].y = cloud->points[i].z;
		cloud_proj->points[i].z = 0;
	}

	cloud_center->points[0].x = mid_y;
	cloud_center->points[0].y = mid_z;
	cloud_center->points[0].z = 0;

	/*


	for (int i = 0; i < cloud->points.size(); i++)
	{
	cloud_tubao_ccw->points[i].x = cloud->points[i].x;
	cloud_tubao_ccw->points[i].y = cloud->points[i].z;
	cloud_tubao_ccw->points[i].z = min_z;
	}

	for (int i = 0; i < cloud->points.size(); i++)
	{
	cloud_tubao_cw->points[i].x = cloud->points[i].x;
	cloud_tubao_cw->points[i].y = cloud->points[i].z;
	cloud_tubao_cw->points[i].z = min_z;
	}



	for (int i = 0; i < cloud->points.size(); i++)
	{
	cloud_proj->points[i].x = cloud->points[i].x;
	cloud_proj->points[i].y = cloud->points[i].z;
	cloud_proj->points[i].z = min_z;
	}

	mpoint *points;
	points = new mpoint[size];
	for (int i = 0; i < cloud->points.size(); i++)
	{
	points[i].x = cloud->points[i].x;
	points[i].y = cloud->points[i].z;
	}
	//**************ccw***************
	int miny_point_id_ccw;
	double *mcos_ccw;
	mcos_ccw = new double[size];
	//point:所需计算的数据计算凸包的存贮数组





	miny_point_id_ccw = get_miny_point_id(points, size);

	get_cos(points, mcos_ccw, miny_point_id_ccw, size);

	sort_points_down(points, mcos_ccw, size);


	get_outpoint_down(points, size, cloud_tubao_ccw);

	//********calc rad**********
	//************ccw**********
	vector<double> rad_ccw;

	for (int i = 0; i < cloud_tubao_ccw->points.size(); i++)
	{
	rad_ccw.push_back(calcAzimuthAngle(mid_x, mid_z, cloud_tubao_ccw->points[i].x, cloud_tubao_ccw->points[i].z));
	}
	cout << "1" << endl;

	//*************cw*************
	int miny_point_id_cw;
	double *mcos_cw;
	mcos_cw = new double[size];
	miny_point_id_cw = get_miny_point_id(points, size);

	get_cos(points, mcos_cw, miny_point_id_cw, size);

	sort_points_up(points, mcos_cw, size);


	get_outpoint_up(points, size, cloud_tubao_cw);

	//********calc rad**********
	//************cw**********
	vector<double> rad_cw;

	for (int i = 0; i < cloud_tubao_cw->points.size(); i++)
	{
	rad_cw.push_back(calcAzimuthAngle(mid_x, mid_z, cloud_tubao_cw->points[i].x, cloud_tubao_cw->points[i].z));
	}


	//vector<double> rad;


	for (int i = 0; i < rad_ccw.size(); i++)
	{
	rad.push_back(rad_ccw[i]);
	}

	for (int i = 0; i < rad_cw.size(); i++)
	{
	rad.push_back(rad_cw[i]);
	}

	cloud_tubao_result->width = (cloud_tubao_cw->points.size() + cloud_tubao_ccw->points.size());
	cloud_tubao_result->height = 1;
	cloud_tubao_result->points.resize(cloud_tubao_result->width * cloud_tubao_result->height);
	cloud_tubao_result->is_dense = true;

	for (int i = 0; i < cloud_tubao_cw->points.size(); i++)
	{
	cloud_tubao_result->points[i].x = cloud_tubao_cw->points[i].x;
	cloud_tubao_result->points[i].y = cloud_tubao_cw->points[i].y;
	cloud_tubao_result->points[i].z = cloud_tubao_cw->points[i].z;
	}

	int j = 0;
	for (int i = cloud_tubao_cw->points.size(); i < (cloud_tubao_cw->points.size() + cloud_tubao_ccw->points.size()); i++)
	{
	cloud_tubao_result->points[i].x = cloud_tubao_ccw->points[j].x;
	cloud_tubao_result->points[i].y = cloud_tubao_ccw->points[j].y;
	cloud_tubao_result->points[i].z = cloud_tubao_ccw->points[j].z;
	j++;
	}

	vector<vector<double> > rad(cloud_tubao_result->points.size(), vector<double>(2));
	int num = 0;
	for (int i = 0; i < cloud_tubao_result->points.size(); i++){
	rad[i][0] = calcAzimuthAngle(mid_x, mid_z, cloud_tubao_result->points[i].x, cloud_tubao_result->points[i].y);
	rad[i][1] = num;
	num++;
	}

	//cout << " before" << endl;
	//for (int i = 0; i < cloud_tubao_result->points.size(); i++)
	//{
	//	for (int j = 0; j < 2; j++)
	//	{
	//		cout << rad[i][j] << " ";
	//	}
	//	cout << endl;
	//}

	sort(rad.begin(), rad.end());

	cout << " after" << endl;
	for (int i = 0; i < cloud_tubao_result->points.size(); i++)
	{
	for (int j = 0; j < 2; j++)
	{
	cout << rad[i][j] << " ";
	}
	cout << endl;
	}


	vector<double>dis;

	for (int i = 0; i < cloud_tubao_result->points.size() - 1; i++)
	{
	dis.push_back(getDistance(cloud_tubao_result->points[i].x, cloud_tubao_result->points[i + 1].x, cloud_tubao_result->points[i].y, cloud_tubao_result->points[i + 1].y));
	}





	for (int i = 0; i < size; i++)
	{
	rad.push_back(calcAzimuthAngle(mid_x, mid_y, cloud_tubao_result->points[i].x, cloud_tubao_result->points[i].y));
	}
	for (int i = 0; i < rad.size(); i++)
	{
	cout << rad[i] << endl;
	}

	//sort(rad.begin(), rad.end());
	//rad.erase(unique(rad.begin(), rad.end()), rad.end());
	//for (int i = 0; i < rad.size(); i++)
	//{
	//	cout << rad[i] << endl;
	//}

	//sort(rad.begin(), rad.end());
	//vector<double>::iterator it, it_flag;
	//vector<int>indx;
	//int ind = 0;
	//for (it = ++rad.begin(); it != rad.end(); ind++)
	//{
	//	it_flag = find(rad.begin(), it, *it);
	//	if (it_flag != it)
	//		it = rad.erase(it);
	//	else{
	//		it++;
	//
	//		indx.push_back(ind);
	//	}
	//
	//}
	//for (int i = 0; i < indx.size(); i++)
	//{
	//	cout << indx[i] << endl;
	//}

	cloud_tubao_final->width = cloud_tubao_result->points.size();
	cloud_tubao_final->height = 1;
	cloud_tubao_final->points.resize(cloud_tubao_final->width * cloud_tubao_final->height);
	cloud_tubao_final->is_dense = true;

	for (int i = 0; i < cloud_tubao_final->points.size(); i++)
	{
	cloud_tubao_final->points[i].x = cloud_tubao_result->points[rad[i][1]].x;
	cloud_tubao_final->points[i].y = cloud_tubao_result->points[rad[i][1]].y;
	cloud_tubao_final->points[i].z = cloud_tubao_result->points[rad[i][1]].z;
	}

	cout << " 2" << endl;

	//cout<<lower_hull.size()<<endl;

	cout << "3 " << endl;


	cout << " 4 " << endl;


	cloud_proj->width = size;
	cloud_proj->height = 1;
	cloud_proj->points.resize(cloud_proj->width * cloud_proj->height);
	cloud_proj->is_dense = true;
	*/

	//计算凸包
	calcConvexHull(y, z, size, cloud_tubao_final);


	//pcl::io::savePCDFileASCII("zhutou3_final.pcd", *cloud_tubao_final);


	//原始点云显示窗口
	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer_ori(new pcl::visualization::PCLVisualizer("原始点云"));
	//viewer_ori->addPointCloud<pcl::PointXYZ>(cloud, "cloud");
	//pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> red(cloud_proj, 255, 0, 0);
	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> red(cloud_center, 255, 0, 0);
	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> green(cloud_tubao_final, 0, 255, 0);
	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> blue(cloud_proj, 0, 0, 255);
	viewer_ori->addPointCloud(cloud_proj, blue, "cloud");
	viewer_ori->addPointCloud(cloud_tubao_final, green, "ccw");
	//viewer_ori->addPointCloud(cloud_tubao_ccw, blue, "res");
	viewer_ori->addPointCloud(cloud_center, red, "center");
	//viewer_ori->addPointCloud(cloud_tubao_result, blue, "result");
	//pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> boundaryoints_color_handler(cloud_boundaries, 0, 255, 0);
	//viewer_ori->addPointCloud(cloud_boundaries, boundaryoints_color_handler, "boundary");
	//viewer_ori->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "cw");
	viewer_ori->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "ccw");
	viewer_ori->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "center");
	viewer_ori->addCoordinateSystem(1);

	/*
	for (int i = 0; i < cloud_tubao_ccw->points.size() - 1; i++)
	{
	//viewer_ori->addText3D(to_string(i), cloud_tubao->points[i], 5.0, 0.5, 0.5, 0.5);
	//viewer_ori->addArrow(cloud_tubao->points[i], cloud_tubao->points[i + 1], 0.5, 0.5, 0.5, to_string(i), false);
	viewer_ori->addLine<pcl::PointXYZ>(cloud_tubao_ccw->points[i], cloud_tubao_ccw->points[i + 1], 0.5, 0.5, 0.5, to_string(i));
	}

	viewer_ori->addLine<pcl::PointXYZ>(cloud_tubao_ccw->points[cloud_tubao_ccw->points.size() - 1], cloud_tubao_ccw->points[0], 0.5, 0.5, 0.5, to_string(cloud_tubao_ccw->points.size() + 1));

	for (int i = 0; i < cloud_tubao_cw->points.size() - 1; i++)
	{
	//viewer_ori->addText3D(to_string(i), cloud_tubao_cw->points[i], 0.5, 0.5, 0.5, 0.5);
	//viewer_ori->addArrow(cloud_tubao->points[i], cloud_tubao->points[i + 1], 0.5, 0.5, 0.5, to_string(i), false);
	viewer_ori->addLine<pcl::PointXYZ>(cloud_tubao_cw->points[i], cloud_tubao_cw->points[i + 1], 0.7, 0.7, 0.7, to_string((i + 100)));
	}

	viewer_ori->addLine<pcl::PointXYZ>(cloud_tubao_cw->points[cloud_tubao_cw->points.size() - 1], cloud_tubao_cw->points[0], 0.7, 0.7, 0.7, to_string(cloud_tubao_cw->points.size() + 1));

	*/

	/*
	for (int i = 0; i < cloud_tubao_final->points.size() - 1; i++)
	{
	//viewer_ori->addText3D(to_string(i), cloud_tubao->points[i], 5.0, 0.5, 0.5, 0.5);
	//viewer_ori->addArrow(cloud_tubao->points[i], cloud_tubao->points[i + 1], 0.5, 0.5, 0.5, to_string(i), false);
	viewer_ori->addLine<pcl::PointXYZ>(cloud_tubao_final->points[i], cloud_tubao_final->points[i + 1], 0.5, 0.5, 0.5, to_string(i));
	}

	viewer_ori->addLine<pcl::PointXYZ>(cloud_tubao_final->points[cloud_tubao_final->points.size() - 1], cloud_tubao_final->points[0], 0.5, 0.5, 0.5, to_string(cloud_tubao_final->points.size() + 1));
	*/
	//showID(cloud_tubao_final, viewer_ori, 1.0, 0.5, 0.5, 0.5);
	drawLine(cloud_tubao_final, viewer_ori, 0.5, 0.5, 0.5);
	cout << "original cloud visulization finish" << endl;


	DWORD end_time = GetTickCount(); //结束计时
	cout << "The run time is:" << (end_time - start_time) << "ms!" << endl; //输出时间


	//viewer->addPointCloudNormals<pcl::PointXYZ, pcl::Normal>(cloud, cloud_normal, 20, 0.03, "normals");
	while (!viewer_ori->wasStopped())
	{
		viewer_ori->spinOnce(100);
		boost::this_thread::sleep(boost::posix_time::microseconds(100000));
	}
}