//Copyright (c) 2016 nbrancati

// The code is free to use for research, provided that the following paper is cited in the works that use the code. 
// Brancati N., De Pietro G., Frucci M., Gallo L., "Dynamic color clustering for human skin detection
// under severe illumination variations", Computer Vision and Image Understanding (submitted)
#include "cv.h"
#include "cvaux.h"
#include "cxmisc.h"
#include "highgui.h"
#include "cxcore.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "color_const.c"


int bins=256;
int tolCr=1;
int tolCb=1;
int useLogTime=1;
float ranges1[] = { 0, 255 };
float ranges2[] = { 0, 255 };
float* ranges[] = { ranges1, ranges2 };
int hist_size[] = {bins, bins};
CvHistogram *histYCr,*histYCb;
int YRCrMin,YRCrMax,YRCbMin,YRCbMax;
char pathf[400]="original_images/";
char *image;
char pathout[400]="result/";
char delimiter[200]=".";
char filename[400]="",immagine[200],output[400]="";

//sort of the histogram
int *sortHist(int *iBins,int *values,int num) {
  int i, tmp, ultimo,itmp;
   int tmpN=num; 
 
   while (tmpN >= 0) 
     { 
         ultimo = -1;
         for (i=0; i<tmpN; i++)
         {
           if (iBins[i]>iBins[i+1]) 
           { 
             tmp = iBins[i]; 
             iBins[i] = iBins[i+1]; 
             iBins[i+1] = tmp;
             ultimo = i;
	     tmp = values[i]; 
             values[i] = values[i+1]; 
             values[i+1] = tmp;
           } 
         }
     tmpN = ultimo;
     }
return iBins;
}


//Computation of Min e max of the histogram
void calcMinMaxHist(int *yValues,int *iBins,int vect[2]){
	int i,j,flag=0,*app,k,y;
	float maxVal=0,percentage=0;
	app=(int*)calloc(bins,sizeof(int));
	for (i=0;i<yValues[0];i++) app[i]=0;
	for (i=1;i<yValues[0];i++){
			maxVal=maxVal+yValues[i];
	}
	i=1;
	if ((int)maxVal!=0){
	while (!flag){
		percentage=percentage+((int)yValues[i]);
		
		if (ceil((percentage/maxVal)*100)>=5) {flag=1;y=yValues[i];}
		i++;
	}
	vect[0]=i-1;
	i=1;flag=0;percentage=0;
	while (!flag){
		percentage=percentage+((int)yValues[i]);
		if (ceil((percentage/maxVal)*100)>=95) {flag=1;y=yValues[i];}
		i++;
	}
	vect[1]=i-1;
	
	k=0;
	for (i=vect[0];i<=vect[1];i++){
		if(iBins[i]!=0){
			app[k]=iBins[i];
			k++;
		}
	}
	app=sortHist(app,iBins,k-1);
	vect[0]=255;vect[1]=0;
	for (i=0;i<k;i++) {
		if (app[i]!=0){vect[0]=app[i];break;}
	}
	for (i=k-1;i>=0;i--) {
		if (app[i]!=0){vect[1]=app[i];break;}
	
	}
	}
	else {vect[0]=255;vect[1]=0;}
}




//Computation of Y0 and Y1 for the computation of the vertices B, C, F and G of the trapezium
void calculateValueMinMaxY(IplImage *image,double val,CvHistogram *hist,int minMax[2],int canal){
int indTol;
int value,i,j,k,**yValue,min=255,max=0,vMin,vMax,*iBins,*app,*app2,*iBins2,indMax=0,indMin=0,**iBinsVal,tol;
uchar spk,l;
double tmpVal=val;
if(canal==1)tol=tolCr;else tol=tolCb;
indTol=(2*(tol+1))-1;
	app=(int*)calloc(bins,sizeof(int));
	iBins=(int*)calloc(bins,sizeof(int));
	app2=(int*)calloc(bins,sizeof(int));
	iBins2=(int*)calloc(bins,sizeof(int));
	for (i=0;i<bins;i++) {app[i]=0;app2[i]=0;iBins2[i]=0;iBins[i]=0;}
	yValue=(int**)calloc(indTol,sizeof(int*));
	iBinsVal=(int**)calloc(indTol,sizeof(int*));
	for (i=0; i<indTol;i++){
		yValue[i]=(int*)calloc(bins,sizeof(int*));
		iBinsVal[i]=(int*)calloc(bins,sizeof(int*));
	}
	for (j=0;j<indTol;j++) for (i=0;i<bins;i++) {yValue[j][i]=0;iBinsVal[j][i]=0;}
	
	for (i=0;i<image->height-1;i++){
		for (j=0;j<image->width-1;j++){

			spk=((uchar *)(image->imageData + i*image->widthStep))[j*image->nChannels + canal];
			if (spk>=tmpVal-tol && spk<=tmpVal+tol) {
				k=((uchar *)(image->imageData + i*image->widthStep))[j*image->nChannels + 0];
                		    int bin_val =0;
				 	bin_val = cvQueryHistValue_2D( hist, k,spk);
				    if (bin_val!=0){
				for (l=0;l<indTol;l++){
						if ((int)(tmpVal-spk+l)==tol) {
							yValue[l][k]=bin_val;iBinsVal[l][k]=k;
						}
					}
        			}
			}
		}
	}
	for (i=0;i<indTol;i++){
		for( k = 0; k < bins; k++ ){app[k]=yValue[i][k];iBins[k]=iBinsVal[i][k];}
			app=sortHist(app,iBins,bins);
			j=1;
			for (k=0;k<bins;k++)if (app[k]!=0){app2[j]=app[k];iBins2[j]=iBins[k];j++;}
			app2[0]=j;
			minMax[0]=255;minMax[1]=0;
			calcMinMaxHist(app2,iBins2,minMax);
			if(minMax[0]!=minMax[1]){	
				if(minMax[0]!=255) {indMin++; if (minMax[0]<min) min=minMax[0];}
				if(minMax[1]!=0) {indMax++; if (minMax[1]>max) max=minMax[1];}
			}
	}
		minMax[0]=min;
		minMax[1]=max;
}

//Computation of histogram of plane1 (Cb and Cr)
int* calculateHist(IplImage *plane1)
{
	int i;
	int *histogram;
	histogram=(int *)calloc(bins,sizeof(int));
	for( i = 0; i < bins; i++ ) histogram[i]=0;
	float ranges1[] = { 0, 255 };
	float* ranges[] = { ranges1};
        IplImage* planes[] = { plane1};
        int hist_size[] = {bins};
        CvHistogram* hist;
        float max_value = 0;

        
        hist = cvCreateHist( 1, hist_size, CV_HIST_ARRAY, ranges, 1 );
        cvCalcHist( planes, hist, 0, 0);
        cvGetMinMaxHistValue( hist, 0, &max_value, 0, 0 );
        for( i = 0; i < bins; i++ ){
		
                float bin_val = cvQueryHistValue_1D( hist, i );
              	histogram[i]=bin_val;
        }
return histogram;
}

CvHistogram *calculateHist2(IplImage *plane1,IplImage *plane2)
{
	int scale=1;
        IplImage* planes[] = { plane1, plane2 };
        IplImage* hist_img=cvCreateImage( cvSize(bins*scale,bins*scale), 8, 3 );
	cvZero( hist_img );
        CvHistogram* hist;
        float max_value = 0;
        int i, j;
        
        hist = cvCreateHist( 2, hist_size, CV_HIST_ARRAY, ranges, 1 );
        cvCalcHist( planes, hist, 0, 0 );
        cvGetMinMaxHistValue( hist, 0, &max_value, 0, 0 );
      
return hist;
}


/*--------------- SKIN SEGMENTATION ---------------*/
int main (int argc, char **argv) {
  CvCapture *capture =0;
	IplImage *source,*frame_rgb, *frame_ycbcr,*bw_ycbcr,*bw_final,*y_plane,*cr_plane,*cb_plane,*grad;
	int i,j,l,k,percentage=50,minMaxCr[2],minMaxCb[2],height,width,perc;
	int *histCr,*histCb;
	
	histCr=(int *)calloc(bins,sizeof(int));
	histCb=(int *)calloc(bins,sizeof(int));
	histYCr = cvCreateHist( 2, hist_size, CV_HIST_ARRAY, ranges, 1 );
	histYCb = cvCreateHist( 2, hist_size, CV_HIST_ARRAY, ranges, 1 );
	for( i = 0; i < bins; i++ ) {histCr[i]=0;histCb[i]=0;}
	float CrMin=133,CrMax=183,CbMin=77,CbMax=128,YMin=0,YMax=255,HCr,HCb,DCr,DCb,D1Cr,D1Cb,dCr,dCb,CbS,radD1Cb;
	double max_valCr,min_valCb,coeffCr,coeffCb;
	DCr=CrMax-CrMin;
	DCb=CbMax-CbMin;

	strcpy(filename,"");
	strcpy(output,"");
	strcpy(immagine,argv[1]);
	strcat(filename,immagine);
	image=strtok(immagine,delimiter);
	printf("working on %s\n",filename);
	strcpy(output,pathout);strcat(output,image);strcat(output,".bmp");
	if( (source = cvLoadImage( filename, -1 )) == 0 ) return -1;		
	if ((int)(source->width)>(int)(source->height)) {height=400;width=604;}
	else {height=604;width=400;}
	frame_rgb = cvCreateImage( cvSize(source->width , source->height) ,
                                     source->depth, source->nChannels );
	frame_ycbcr = cvCreateImage( cvSize(frame_rgb->width , frame_rgb->height),frame_rgb->depth, frame_rgb->nChannels );
	grad = cvCreateImage( cvSize(frame_rgb->width , frame_rgb->height),8, 1);
	bw_final = cvCreateImage( cvSize(frame_rgb->width , frame_rgb->height),8, 1);	
	y_plane = cvCreateImage( cvSize(frame_rgb->width , frame_rgb->height), 8, 1 );
        cr_plane = cvCreateImage( cvSize(frame_rgb->width , frame_rgb->height), 8, 1 );
        cb_plane = cvCreateImage( cvSize(frame_rgb->width , frame_rgb->height), 8, 1 );
	cvZero(y_plane);
	cvCopy(source,frame_rgb,0);
	perc=frame_rgb->width*frame_rgb->height*0.1/100;

	cvCvtColor(frame_rgb, frame_ycbcr, CV_BGR2YCrCb);
	cvSplit( frame_ycbcr, y_plane, cr_plane, cb_plane, 0 );	

	histCb=calculateHist(cb_plane);
	histCr=calculateHist(cr_plane);

	
	
	max_valCr=0;
	minMaxCr[0]=255;minMaxCr[1]=0;
	minMaxCb[0]=255;minMaxCb[1]=0;

	//Computation of Crmax
	for (i=bins-1;i>=0;i--){
		if (histCr[i]!=0 && histCr[i]>perc) {
				max_valCr=i;
				break;
		}
	}
	
	//Computation of Cbmin
	min_valCb=0;
	for (i=0;i<bins;i++){
		if (histCb[i]!=0 && histCb[i]>perc) {
			min_valCb=i;
			break;
		}
	}



	histYCb=calculateHist2(y_plane,cb_plane);
	histYCr=calculateHist2(y_plane,cr_plane);

//Computation of B and C
	if (max_valCr!=-1) {
		if (max_valCr>CrMax)max_valCr=CrMax;
		calculateValueMinMaxY(frame_ycbcr,max_valCr,histYCr,minMaxCr,1);
		if (max_valCr<CrMax)CrMax=max_valCr;
		
	}	


	//Computation of F and G
	if (min_valCb!=-1){
		if (min_valCb<CbMin)min_valCb=CbMin;
		calculateValueMinMaxY(frame_ycbcr,min_valCb,histYCb,minMaxCb,2);
		if (min_valCb>CbMin)CbMin=min_valCb;
		
	}	


        //
	YRCrMin=50,YRCrMax=110,YRCbMin=140,YRCbMax=200;
	if (max_valCr!=-1){
		YRCrMin=minMaxCr[0];
		YRCrMax=minMaxCr[1];
	}
	if (max_valCr!=-1){
		YRCbMin=minMaxCb[0];
		YRCbMax=minMaxCb[1];
	}
	cvZero(bw_final);
	int Y,Cr,Cb,valueY;


	for (i=0;i<frame_rgb->height-1;i++){
		for (j=0;j<frame_rgb->width-1;j++){
			uchar Y=((uchar *)(frame_ycbcr->imageData + i*frame_ycbcr->widthStep))[j*frame_ycbcr->nChannels + 0];
			uchar Cr=((uchar *)(frame_ycbcr->imageData + i*frame_ycbcr->widthStep))[j*frame_ycbcr->nChannels + 1];
			uchar Cb= ((uchar *)(frame_ycbcr->imageData + i*frame_ycbcr->widthStep))[j*frame_ycbcr->nChannels + 2];
			HCr=0;HCb=0;CbS=0;
			//Pixels in the trapezium
			if ((Cr<=CrMax && Cr>=CrMin) && (Cb<=CbMax && Cb>=CbMin) && Cr>Cb && Y>=YMin && Y<=YMax)
				{//Calculate HCr
				if (Y>=YMin && Y<YRCrMin) {
					coeffCr=(CrMax-CrMin)/(YRCrMin-YMin);
					coeffCr = (round(coeffCr * 100))/100;
					HCr=round(Y*coeffCr+CrMin-YMin*coeffCr);
				}
				else if (Y>=YRCrMin && Y<YRCrMax) HCr=round(CrMax);
				else if (Y>=YRCrMax && Y<=YMax) {
					coeffCr=(CrMin-CrMax)/(YMax-YRCrMax);
					coeffCr = (round(coeffCr * 100))/100;
					HCr=round(Y*coeffCr+CrMax-YRCrMax*coeffCr);
				}
				//Calculate HCb
				if (Y>=YMin && Y<YRCbMin) {
					coeffCb=(CbMin-CbMax)/(YRCbMin-YMin);
					coeffCb = (round(coeffCb * 100))/100;
					HCb=round(Y*coeffCb+CbMax-YMin*coeffCb);
				}
				else if (Y>=YRCbMin && Y<YRCbMax) HCb=round(CbMin);
				else if (Y>=YRCbMax && Y<=YMax) {
					coeffCb=(CbMax-CbMin)/(YMax-YRCbMax);	
					coeffCb = (round(coeffCb * 100))/100;
					HCb=round(Y*coeffCb+CbMin-YRCbMax*coeffCb);
				}

				if (HCr>CrMax)HCr=CrMax;
				if (HCb<CbMin)HCb=CbMin;
				//condition C.2
				dCr=Cr-CrMin;
				D1Cr=HCr-CrMin;
				D1Cb=CbMax-HCb;
				dCb=dCr*D1Cb/D1Cr;
				CbS=round(CbMax-dCb);
				radD1Cb=round(D1Cb/3);
				//Conditions C.0, C.1 and C.2
				if (Cr<=HCr && Cb<=CbS+radD1Cb && Cb>=CbS-radD1Cb && Cb>=HCb)
					cvSet2D(bw_final,i,j,cvScalarAll(255));
		}}}

IplConvKernel *element;
element = cvCreateStructuringElementEx(7, 7, 3, 3, CV_SHAPE_CROSS,0);
cvDilate(bw_final,bw_final,element,1);
cvErode(bw_final,bw_final,element,1);
cvNot(bw_final,bw_final);

	cvSaveImage(output,bw_final);

		
	cvReleaseImage( &frame_rgb );
	cvReleaseImage( &frame_ycbcr );
	cvReleaseImage( &bw_final );

}
