#include "ImageIOpfm.h"
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <iostream>
#include <dirent.h>
#include <fstream>
#include <sstream>
#include <string>

#define PI 3.1415926
#define CIRCLR_ERR 3
#define R 1.5

//amy

cv::Mat K( cv::Matx33f(268.5119, 0, 320,
						 0, 268.5119, 240,
						 0, 0, 1) );

int k = 0;
int m = 0;


void GetFileNames(string path,vector<string>& filenames)
{
    /*DIR *pDir;
    struct dirent* ptr;
    if(!(pDir = opendir(path.c_str())))
        return;
    while((ptr = readdir(pDir))!=0) {
        if (strcmp(ptr->d_name, ".") != 0 && strcmp(ptr->d_name, "..") != 0)
        {
            cout << ptr->d_off << " " << ptr->d_name << endl;
            filenames.push_back(path + "/" + ptr->d_name);
        }
    }
    closedir(pDir);*/
    struct dirent** namelist;
    int n;
    n = scandir(path.c_str(), &namelist, 0, alphasort);
    if(n < 0)
        cout << "scandir return " << n << endl;
    else
    {
        int index = 2;
        while(index < n)
        {
            //cout << namelist[index]->d_name << endl;
            filenames.push_back(path + "/" + namelist[index]->d_name);
            free(namelist[index]);
            index++;
        }
    }
    free(namelist);
}

int otsuThreshold(IplImage* img)
{
	
	int T = 0;
	int height = img->height;
	int width  = img->width;
	int step      = img->widthStep;
	int channels  = img->nChannels;
	uchar* data  = (uchar*)img->imageData;
	double gSum0;
	double gSum1;
	double N0 = 0;
	double N1 = 0;
	double u0 = 0;
	double u1 = 0;
	double w0 = 0;
	double w1 = 0;
	double u = 0;
	double tempg = -1;
	double g = -1;
	double Histogram[256]={0};// = new double[256];
	double N = width*height;
	for(int i=0;i<height;i++)
	{
		for(int j=0;j<width;j++)
		{
			double temp =data[i*step + j * 3] * 0.114 + data[i*step + j * 3+1] * 0.587 + data[i*step + j * 3+2] * 0.299;
			temp = temp<0? 0:temp;
			temp = temp>255? 255:temp;
			Histogram[(int)temp]++;
		} 
	}
	
	for (int i = 0;i<256;i++)
	{
		gSum0 = 0;
		gSum1 = 0;
		N0 += Histogram[i];			
		N1 = N-N0;
		if(0==N1)break;
		w0 = N0/N;
		w1 = 1-w0;
		for (int j = 0;j<=i;j++)
		{
			gSum0 += j*Histogram[j];
		}
		u0 = gSum0/N0;
		for(int k = i+1;k<256;k++)
		{
			gSum1 += k*Histogram[k];
		}
		u1 = gSum1/N1;
		//u = w0*u0 + w1*u1;
		g = w0*w1*(u0-u1)*(u0-u1);
		if (tempg<g)
		{
			tempg = g;
			T = i;
		}
	}
	return T; 
}

double angle( Point pt1, Point pt2, Point pt0 )
{
    double dx1 = pt1.x - pt0.x;
    double dy1 = pt1.y - pt0.y;
    double dx2 = pt2.x - pt0.x;
    double dy2 = pt2.y - pt0.y;
    return (dx1*dx2 + dy1*dy2)/sqrt((dx1*dx1 + dy1*dy1)*(dx2*dx2 + dy2*dy2) + 1e-10);
}

void drawSquares( Mat& image, const vector<vector<Point> >& squares )
{
    for( size_t i = 0; i < squares.size(); i++ )
    {
        const Point* p = &squares[i][0];

        int n = (int)squares[i].size();
        //dont detect the border
        polylines(image, &p, &n, 1, true, Scalar( 255, 0, 0 ), 2, LINE_AA);
    }

}

void RankPoint( vector<Point> &p, int num = 4 )
{
	vector<double> cosr;
	vector<Point> result;
	double ab, abs_a, abs_b;   //vector b is (1, 1) 
	int temp2 = p[0].x - p[0].y; 
	int right_up = 0;

    Point center;
    for(int i = 0; i < num; i++)
    {
        center.x = center.x + p[i].x;
        center.y = center.y + p[i].y;
    }
    center.x = center.x / num;
    center.y = center.y / num;

	//find out the right_up point as first point
	for( int i = 1; i < num; i++)
	{
		if( p[i].x - p[i].y >= temp2 && p[i].y < center.y )
		{
			temp2 = p[i].x - p[i].y;
			right_up = i;
		}
	}
	result.push_back(p[right_up]);
	p.erase( p.begin() + right_up );

	//compute the angel between vector a and (1,1)
	for( int i = 0; i < num - 1; i++ )
	{
		ab = p[i].x - result[0].x + p[i].y - result[0].y;
		abs_a = sqrt((p[i].x - result[0].x) * (p[i].x - result[0].x) + (p[i].y - result[0].y) * (p[i].y - result[0].y));
		abs_b = sqrt(2);
		cosr.push_back(ab / abs_a /abs_b);
	}

	//sort the points' angels from big to small
	int index[num-1];
	int temp;
	for(int i = 0; i < num - 1; i++)
		index[i] = i;

	for(int i = 0; i < num -1; i++)
	{
		for(int j = i; j < num - 1; j++)
		{
			if(cosr[index[i]] < cosr[index[j]])
			{
				temp = index[i];
				index[i] = index[j];
				index[j] = temp;
			}
		}
	}

	//the points are clockwise now
	for(int i = 0; i < num - 1; i++ )
		result.push_back(p[index[i]]);
	
	swap(p,result);
}

int main(int argc, char ** argv){
    vector<string> filename;
    GetFileNames("../../dataset/pics08111347", filename);
    //GetFileNames("../../dataset/depth", filename);

    double time_start, time_duration;
    for(int i = 0; i < filename.size(); i++)
    {
    time_start = clock();
    //ReadFilePFM(I, filename[i]);
    //Mat srcImage = imread("../1.png", IMREAD_COLOR);
    Mat srcImage = imread(filename[i], IMREAD_COLOR);
    imshow("result", srcImage);
    Mat srcImage3, srcImage4;
    srcImage.copyTo(srcImage3);
    srcImage.copyTo(srcImage4);
//    IplImage srcImage2 = IplImage(srcImage);
//    int adaptThresh = otsuThreshold(&srcImage2);
    Mat gray, bImage;
    cvtColor(srcImage, gray, COLOR_BGR2GRAY);

    threshold(gray, bImage, 170, 255, CV_THRESH_BINARY);
    Mat ele = getStructuringElement(MORPH_RECT, Size(3,3));
    //erode(bImage, bImage, ele);
    dilate(bImage, bImage, ele);

    Mat fhsv, special;
    int minh = 0, maxh = 10;
    int mins = 58, maxs = 255;
    int minv = 0, maxv = 100;
    cvtColor(srcImage3, fhsv,COLOR_BGR2HSV); 
    inRange(fhsv, Scalar(minh,mins,minv), Scalar(maxh,maxs,maxv), special);
    dilate(special, special, ele);
    GaussianBlur(special, special, Size(9, 9), 2);
    vector<Vec3f> circles3;
    HoughCircles(special, circles3, HOUGH_GRADIENT, 1, 10, 80, 70, 0, 600);

    GaussianBlur(gray, gray, Size(9, 9), 2);
    vector<Vec3f> circles;
    HoughCircles(gray, circles, HOUGH_GRADIENT, 1, 10, 80, 70, 0, 600);

    vector< vector<Point> > contours;
	vector<Vec4i> hierarchy;
    findContours(bImage, contours, hierarchy, CV_RETR_TREE, CV_CHAIN_APPROX_NONE);
    vector< vector<Point> > squares;
    vector<Point> approx;
    for(int i = 0; i < contours.size(); i++)
    {
        approxPolyDP(contours[i], approx, arcLength(Mat(contours[i]), true)*0.02f, true);
        if( approx.size() == 4 && fabs(contourArea(Mat(approx))) > 200 && isContourConvex(Mat(approx)) )
        {
            double maxCosine = 0;

            for( int j = 2; j < 5; j++ )
            {
                    // find the maximum cosine of the angle between joint edges
                double cosine = fabs(angle(approx[j%4], approx[j-2], approx[j-1]));
                maxCosine = MAX(maxCosine, cosine);
            }

                // if cosines of all angles are small
                // (all angles are ~90 degree) then write quandrange
                // vertices to resultant sequence
            if( maxCosine < 0.7 ) // 0.3
            //cv::RotatedRect &&rect = cv::minAreaRect(contours[i]);
            squares.push_back(approx);
        }
    }

    vector<Mat> num;
    vector<Mat> num_resize;
    for(int i = 0; i < squares.size(); i++)
    {
        //cv::RotatedRect rect = cv::minAreaRect(squares[i]);
        Rect rect = boundingRect(squares[i]);
        //cout << "before: " << squares[i] << endl;
        RankPoint(squares[i]);
        //cout << "after: " << squares[i] << endl;
        //cv::rectangle(srcImage3, rect, (255, 255, 255), 1);
        Mat imageROI;
        imageROI = srcImage3(rect);
        num.push_back(imageROI);
        /*cv::Point2f pts_src[] = { 
				squares[i][3],
				squares[i][0],
				squares[i][1],
				squares[i][2]
			};
        cv::Point2f pts_dst[] = { 
				cv::Point(rect.tl().x, rect.tl().y),
				cv::Point(rect.tl().x + rect.width, rect.tl().y),
				cv::Point(rect.tl().x + rect.width, rect.tl().y + rect.height) ,
				cv::Point(rect.tl().x, rect.tl().y + rect.height)
			};*/
        cv::Point2f pts_src[] = { 
				squares[i][3] - rect.tl(),
				squares[i][0] - rect.tl(),
				squares[i][1] - rect.tl(),
				squares[i][2] - rect.tl()
			};
        cv::Point2f pts_dst[] = { 
				cv::Point(0, 0),
				cv::Point(rect.width, 0),
				cv::Point(rect.width, rect.height) ,
				cv::Point(0, rect.height)
			};
        cv::Mat M = cv::getPerspectiveTransform(pts_src, pts_dst);
        cv::Mat warp;
        cv::warpPerspective(imageROI, warp, M, imageROI.size());
        num_resize.push_back(warp);
        //cv::warpPerspective(imageROI, warp, M, srcImage3.size(), cv::INTER_LINEAR + cv::WARP_INVERSE_MAP, cv::BORDER_REPLICATE);
    }

    for(int i = 0; i < num_resize.size(); i++)
    {
        k++;
        stringstream ss;
        ss << (double)clock() / CLOCKS_PER_SEC*1000;
        string img_name = "../data/resize/" + ss.str() + ".png";
        imwrite(img_name, num_resize[i]);
    }

    for(int i = 0; i < num.size(); i++)
    {
        m++;
        stringstream ss2;
        //ss2 << k;
        ss2 << (double) clock() / CLOCKS_PER_SEC * 1000;
        string img_name2 = "../data/original/" + ss2.str() + ".png";
        imwrite(img_name2, num[i]);
    }


    /*for(int i = 0; i < num.size(); i++)
    {
        vector< vector<Point> > num_contours;
        vector<Vec4i> num_hierarchy;
        Mat num_gray, bnum;
        cvtColor(num[i], num_gray, COLOR_RGB2GRAY);
        threshold(num_gray, bnum, 160, 255, CV_THRESH_BINARY);
        findContours(bnum, num_contours, hierarchy, CV_RETR_LIST, CV_CHAIN_APPROX_NONE);
        cout << "num contours: " << num_contours.size() << endl;
        vector< vector<Point> >::iterator itc = num_contours.begin();
        while(itc != num_contours.end())
        {
        Rect num_rect = boundingRect(*itc);   //Get the bound rect of the contour
        if(itc->size() < 50)
        {
            itc++;
            continue;
        }
        Mat z;
        bnum.copyTo(z);
        Mat num_imageROI;
        num_imageROI = z(num_rect);
        itc++;
        imshow("test", num_imageROI);
        if(waitKey(0) == 32)
            continue;
        //num.push_back(num_imageROI);
        }
    }*/

    drawSquares(srcImage, squares);
    drawSquares(srcImage3, squares);
    
    vector<Point3f> Pcs;
    Point3f Pc_temp;
    for(int i = 0; i < circles.size(); i++)
    {
        Point center(cvRound(circles[i][0]), cvRound(circles[i][1]));
		int radius = cvRound(circles[i][2]);

        Pc_temp.z = K.at<float>(0,0) * R / radius;
        Pc_temp.x = (center.x - K.at<float>(0,2)) * Pc_temp.z / K.at<float>(0,0);
	    Pc_temp.y = (center.y - K.at<float>(1,2)) * Pc_temp.z / K.at<float>(1,1);
        Pcs.push_back(Pc_temp);

        cout << Pc_temp << endl;

        circle(srcImage, center, 3, Scalar(255, 0, 0), -1, 8, 0);
        circle(srcImage, center, radius, Scalar(0, 255, 0), 3, 8, 0);
    }

    vector<Point3f> Pcs3;
    Point3f Pc_temp3;
    for(int i = 0; i < circles3.size(); i++)
    {
        Point center3(cvRound(circles3[i][0]), cvRound(circles3[i][1]));
		int radius3 = cvRound(circles3[i][2]);

        Pc_temp3.z = K.at<float>(0,0) * R / radius3;
        Pc_temp3.x = (center3.x - K.at<float>(0,2)) * Pc_temp3.z / K.at<float>(0,0);
	    Pc_temp3.y = (center3.y - K.at<float>(1,2)) * Pc_temp3.z / K.at<float>(1,1);
        Pcs3.push_back(Pc_temp3);

        cout << Pc_temp3 << endl;

        circle(srcImage3, center3, 3, Scalar(255, 0, 0), -1, 8, 0);
        circle(srcImage3, center3, radius3, Scalar(0, 255, 0), 3, 8, 0);
    }
   
    //imshow("result", srcImage);
    imshow("bImage", bImage);
    imshow("special", special);
    imshow("result3", srcImage3);

    Pcs.clear();
    Pcs3.clear();

    time_duration = (clock() - time_start)/CLOCKS_PER_SEC*1000;
    cout << "time cost: " << time_duration << endl;
    
    //if(waitKey(0) == 32)
    //    continue;
    waitKey(1);
    }
    return 0;
}
