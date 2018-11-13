#ifndef CONVEXHULL_H
#define CONVEXHULL_H

#include <iostream>
#include <cstdio>
#include <algorithm>
#include <list>

#include <pcl/io/io.h>
#include <pcl/io/pcd_io.h>		//PCL��PCD��ʽ�ļ����������ͷ�ļ�
#include <pcl/io/obj_io.h>
#include <pcl/PolygonMesh.h>
#include <pcl/point_cloud.h>
#include <pcl/io/vtk_lib_io.h>//loadPolygonFileOBJ����ͷ�ļ���
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/sample_consensus/method_types.h>   //����������Ʒ���ͷ�ļ�
#include <pcl/sample_consensus/model_types.h>   //ģ�Ͷ���ͷ�ļ�
#include <pcl/segmentation/sac_segmentation.h>   //���ڲ���һ���Էָ�����ͷ�ļ�
#include <pcl/common/impl/io.hpp>
#include <pcl/point_types.h>
#include <pcl/ModelCoefficients.h>
#include <pcl/filters/project_inliers.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/filters/voxel_grid.h>
#include <boost/thread/thread.hpp>
#include <pcl/features/normal_3d_omp.h>
#include <pcl/features/normal_3d.h>
#include <pcl/features/boundary.h>			//�߽���ȡͷ�ļ�	

#include "myPoint.h"

using namespace std;
using namespace pcl;





class getConvexhull{


public:

	getConvexhull();
	virtual ~getConvexhull();


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
	int get_miny_point_id(mpoint *points, int size);

	//������ֵ
	//************************************
	// Method:    get_cos
	// FullName:  get_cos
	// Access:    public 
	// Returns:   void
	// Qualifier:
	// Parameter: mpoint * points	���е�����
	// Parameter: double * mcos		���е�����ֵ����
	// Parameter: int id			��Сyֵ����
	// Parameter: int size			�����С
	//************************************
	void get_cos(mpoint *points, double *mcos, int id, int size);

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
	double calcAzimuthAngle(double x0, double y0, double x, double y);

	//��������ֵ�Ӵ�С����
	//************************************
	// Method:    sort_points_down
	// FullName:  sort_points_down
	// Access:    public 
	// Returns:   void
	// Qualifier:
	// Parameter: mpoint * points		���е�����
	// Parameter: double * mcos			���е�����ֵ����
	// Parameter: int size				�����С
	//************************************
	void sort_points_down(mpoint *points, double *mcos, int size);

	//������ʱ�뷽��ĵ�
	//************************************
	// Method:    ccw
	// FullName:  ccw
	// Access:    public 
	// Returns:   int
	// Qualifier:
	// Parameter: mpoint a		��ʼ��
	// Parameter: mpoint b		�м��
	// Parameter: mpoint c		ĩβ��
	//************************************
	int ccw(mpoint a, mpoint b, mpoint c);

	//��������ֵ��С��������
	//************************************
	// Method:    sort_points_up
	// FullName:  sort_points_up
	// Access:    public 
	// Returns:   void
	// Qualifier:
	// Parameter: mpoint * points		���е�����
	// Parameter: double * mcos			���е�����ֵ����
	// Parameter: int size				�����С
	//************************************
	void sort_points_up(mpoint *points, double *mcos, int size);

	//����˳ʱ�뷽��ĵ�
	//************************************
	// Method:    cw
	// FullName:  cw
	// Access:    public 
	// Returns:   int
	// Qualifier:
	// Parameter: mpoint a		��ʼ��
	// Parameter: mpoint b		�м��
	// Parameter: mpoint c		ĩβ��
	//************************************
	int cw(mpoint a, mpoint b, mpoint c);

	//�õ�ccw�Ľ����
	//************************************
	// Method:    get_outpoint_down
	// FullName:  get_outpoint_down
	// Access:    public 
	// Returns:   void
	// Qualifier:
	// Parameter: mpoint * points			���е�����
	// Parameter: int size					�����С
	// Parameter: pcl::PointCloud<pcl::PointXYZ>::Ptr & cloud_tubao		�������
	//************************************
	void get_outpoint_down(mpoint *points, int size, pcl::PointCloud<pcl::PointXYZ>::Ptr& cloud_tubao);

	//�õ�cw�Ľ����
	//************************************
	// Method:    get_outpoint_up
	// FullName:  get_outpoint_up
	// Access:    public 
	// Returns:   void
	// Qualifier:	
	// Parameter: mpoint * points		���е�����
	// Parameter: int size				�����С
	// Parameter: pcl::PointCloud<pcl::PointXYZ>::Ptr & cloud_tubao		�������
	//************************************
	void get_outpoint_up(mpoint *points, int size, pcl::PointCloud<pcl::PointXYZ>::Ptr& cloud_tubao);


	//��������Сֵ
	//************************************
	// Method:    min_array
	// FullName:  min_array
	// Access:    public 
	// Returns:   double
	// Qualifier: 
	// Parameter: double * array	��Ҫ�������
	// Parameter: int n				Ԫ�ظ���
	//************************************
	double
	min_array(double *array, int n);

	//���������ֵ
	//************************************
	// Method:    max_array
	// FullName:  max_array
	// Access:    public 
	// Returns:   double
	// Qualifier:
	// Parameter: double * array	��Ҫ�������
	// Parameter: int n				Ԫ�ظ���
	//************************************
	double
	max_array(double *array, int n);

	//���������
	//************************************
	// Method:    getDistance
	// FullName:  getDistance
	// Access:    public 
	// Returns:   double
	// Qualifier:
	// Parameter: double x1		��һ��x
	// Parameter: double y1		��һ��y
	// Parameter: double x2		�ڶ���x
	// Parameter: double y2		�ڶ���y
	//************************************
	double getDistance(double x1, double y1, double x2, double y2);

	//��ʾ���
	//************************************
	// Method:    showID
	// FullName:  showID
	// Access:    public 
	// Returns:   void
	// Qualifier:
	// Parameter: pcl::PointCloud<pcl::PointXYZ>::Ptr & cloud			��ʾ����
	// Parameter: boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer	��ʾ����
	// Parameter: double scale		�ֺ�
	// Parameter: double r			��ɫr
	// Parameter: double g			��ɫg
	// Parameter: double b			��ɫb
	//************************************
	void showID(pcl::PointCloud<pcl::PointXYZ>::Ptr& cloud,
	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer,
	double scale,
	double r,
	double g,
	double b);

	//��˳��������
	//************************************
	// Method:    drawLine
	// FullName:  drawLine
	// Access:    public 
	// Returns:   void
	// Qualifier:
	// Parameter: pcl::PointCloud<pcl::PointXYZ>::Ptr & cloud		��ʾ����
	// Parameter: boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer		��ʾ����
	// Parameter: double r			��ɫr
	// Parameter: double g			��ɫg
	// Parameter: double b			��ɫb
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
	// Parameter: double * x	ͶӰ����x��
	// Parameter: double * y	ͶӰ����y��
	// Parameter: int size		���ƴ�С
	// Parameter: pcl::PointCloud<pcl::PointXYZ>::Ptr & cloud_final	���ս������
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
	// Parameter: string file	�ļ�����(·��)
	//************************************
	void ch_operation(string file);

};
#endif	//CONVEXHULL_H
