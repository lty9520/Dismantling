#ifndef CONVEXHULL_H
#define CONVEXHULL_H

#include <iostream>
#include <cstdio>
#include <algorithm>
#include <list>

#include <pcl/io/io.h>
#include <pcl/io/pcd_io.h>		//PCL的PCD格式文件的输入输出头文件
#include <pcl/io/obj_io.h>
#include <pcl/PolygonMesh.h>
#include <pcl/point_cloud.h>
#include <pcl/io/vtk_lib_io.h>//loadPolygonFileOBJ所属头文件；
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/sample_consensus/method_types.h>   //随机参数估计方法头文件
#include <pcl/sample_consensus/model_types.h>   //模型定义头文件
#include <pcl/segmentation/sac_segmentation.h>   //基于采样一致性分割的类的头文件
#include <pcl/common/impl/io.hpp>
#include <pcl/point_types.h>
#include <pcl/ModelCoefficients.h>
#include <pcl/filters/project_inliers.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/filters/voxel_grid.h>
#include <boost/thread/thread.hpp>
#include <pcl/features/normal_3d_omp.h>
#include <pcl/features/normal_3d.h>
#include <pcl/features/boundary.h>			//边界提取头文件	

#include "myPoint.h"

using namespace std;
using namespace pcl;





class getConvexhull{


public:

	getConvexhull();
	virtual ~getConvexhull();


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
	int get_miny_point_id(mpoint *points, int size);

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
	void get_cos(mpoint *points, double *mcos, int id, int size);

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
	double calcAzimuthAngle(double x0, double y0, double x, double y);

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
	void sort_points_down(mpoint *points, double *mcos, int size);

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
	int ccw(mpoint a, mpoint b, mpoint c);

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
	void sort_points_up(mpoint *points, double *mcos, int size);

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
	int cw(mpoint a, mpoint b, mpoint c);

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
	void get_outpoint_down(mpoint *points, int size, pcl::PointCloud<pcl::PointXYZ>::Ptr& cloud_tubao);

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
	void get_outpoint_up(mpoint *points, int size, pcl::PointCloud<pcl::PointXYZ>::Ptr& cloud_tubao);


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
	min_array(double *array, int n);

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
	max_array(double *array, int n);

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
	double getDistance(double x1, double y1, double x2, double y2);

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
	void showID(pcl::PointCloud<pcl::PointXYZ>::Ptr& cloud,
	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer,
	double scale,
	double r,
	double g,
	double b);

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
	void drawLine(pcl::PointCloud<pcl::PointXYZ>::Ptr& cloud,
	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer,
	double r,
	double g,
	double b);

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
	void calcConvexHull(double *x,
	double *y,
	int size,
	pcl::PointCloud<pcl::PointXYZ>::Ptr& cloud_final);

	//************************************
	// Method:    ch_operation
	// FullName:  ch_operation
	// Access:    public 
	// Returns:   void
	// Qualifier:
	// Parameter: string file	文件名称(路径)
	//************************************
	void ch_operation(string file);

};
#endif	//CONVEXHULL_H
