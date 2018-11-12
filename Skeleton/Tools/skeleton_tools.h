#ifndef SKELETON_H
#define SKELETON_H

#include <opencv2/imgproc/imgproc.hpp>  
#include <opencv2/objdetect/objdetect.hpp>  
#include <opencv2/highgui/highgui.hpp>  

#include <list>
#include <set>
#include <iostream> 
#include <algorithm>  

using namespace std;
using namespace cv;

class skeletonThinning
{

private:
	


public:
	skeletonThinning();
	virtual ~skeletonThinning();
	//获取A0~A5
	//************************************
	// Method:    GetAi
	// FullName:  GetAi
	// Access:    public 
	// Returns:   set<int>
	// Qualifier:
	// Parameter: int a[]
	// Parameter: int length
	//************************************
	set<int> GetAi(int a[], int length);


	//迭代腐蚀
	//************************************
	// Method:    erodephase
	// FullName:  erodephase
	// Access:    public 
	// Returns:   bool
	// Qualifier:
	// Parameter: list<cv::Point> & border
	// Parameter: cv::Mat & Input
	// Parameter: int neighbour[][3]
	// Parameter: const set<int> & A
	//************************************
	bool erodephase(list<cv::Point> &border, cv::Mat&Input, int neighbour[][3], const set<int>& A);


	//找边界 
	//************************************
	// Method:    findborder
	// FullName:  findborder
	// Access:    public 
	// Returns:   void
	// Qualifier:
	// Parameter: list<cv::Point2i> & border
	// Parameter: const cv::Mat & Input
	//************************************
	void findborder(list<cv::Point2i>& border, const cv::Mat&Input);


	//最后一步，得到骨架
	//************************************
	// Method:    finalerode
	// FullName:  finalerode
	// Access:    public 
	// Returns:   void
	// Qualifier:
	// Parameter: cv::Mat & Input
	// Parameter: int neighbour[][3]
	// Parameter: const set<int> & A
	//************************************
	void finalerode(cv::Mat&Input, int neighbour[][3], const set<int>& A);


	//************************************
	// Method:    thin
	// FullName:  thin
	// Access:    public 
	// Returns:   void
	// Qualifier: //Input是二值图像
	// Parameter: cv::Mat & Input
	//************************************
	void thin(cv::Mat &Input);


	/**
	* Perform one thinning iteration.
	* Normally you wouldn't call this function directly from your code.
	*
	* Parameters:
	* 		im    Binary image with range = [0,1]
	* 		iter  0=even, 1=odd
	*/
	void thinningIteration(cv::Mat& img, int iter);


	/**
	* Function for thinning the given binary image
	*
	* Parameters:
	* 		src  The source image, binary with range = [0,255]
	* 		dst  The destination image
	*/
	void thinning(const cv::Mat& src, cv::Mat& dst);

	//************************************
	// Method:    zhang_thinimage_improve
	// FullName:  zhang_thinimage_improve
	// Access:    public 
	// Returns:   void
	// Qualifier: //单通道、二值化后的图像
	// Parameter: Mat & srcimage
	//************************************
	void zhang_thinimage_improve(Mat &srcimage);

	

};


#endif //SKELETON_H
