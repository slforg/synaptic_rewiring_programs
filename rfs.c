#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

#include<X11/Xlib.h>
#include<X11/Xutil.h>
#include<unistd.h>

#define XSIZE 24
#define YSIZE 24
#define COOR_MAX 256
#define WINDOW_X 24
#define WINDOW_Y 24
#define DOT_SIZE 12
#define WINDOW_WIDTH XSIZE*DOT_SIZE
#define WINDOW_HEIGHT YSIZE*DOT_SIZE
#define WINDOW_F 10
#define istm0 27
#define istm00 27

//Color Return Function ex GetColor(dis,"yellow");
unsigned long GetColor(Display* dis, char* color_name){
  Colormap cmap;
  XColor near_color,true_color;

  cmap=DefaultColormap(dis,0);
  XAllocNamedColor(dis,cmap,color_name,&near_color,&true_color);
  return (near_color.pixel);
}

//Dot Color
void DotRGB(Display *dis,Window win,GC gc,int x,int y,unsigned char r,unsigned char g,unsigned char b){  XSetForeground(dis,gc,(r<<16)|(g<<8)|b);
  XFillRectangle(dis,win,gc,x,y,DOT_SIZE,DOT_SIZE);
}

main(){
  int x,y,x0,y0,dx,dy,xd,yd,i,j,i0,j0,k0;
  int sx,sy; //select x&y point
  int dir_i;
  char name[10];
  int color_Set[256][256];
  float dir[XSIZE][YSIZE][istm0],rmx,rmn,nmlz,rnm,rr;

  Display *dis;
  Window rw,w;
  GC gc;
  XSetWindowAttributes att;

  typedef struct{
    unsigned int r,g,b;
  }COLOR;

  COLOR rgb[256];

  FILE *fdir;

  printf("input file name : ");
  scanf("%s",name);
  getchar();

  if((fdir = fopen(name,"rb"))==NULL){
    printf("Can not open output file.\n");
    exit(-1);
  }
  for(k0=0;k0<istm0;k0++){
    for(j0=0;j0<YSIZE;j0++){
      for(i0=0;i0<YSIZE;i0++){
	fscanf(fdir,"%f\n" ,&dir[i0][j0][k0]);
      }
    }
  }

  dis = XOpenDisplay(NULL);
  rw =XDefaultRootWindow(dis);
 
  w=XCreateSimpleWindow(dis,
			rw,
			WINDOW_X,
			WINDOW_Y,
			WINDOW_WIDTH,
			WINDOW_HEIGHT,
			WINDOW_F,
			WhitePixel(dis,0),
			BlackPixel(dis,0));
  XStoreName(dis,w,"Direction_Map"); //window name
  att.backing_store=Always;
  XChangeWindowAttributes(dis,w,CWBackingStore,&att);
  XMapWindow(dis,w); //Window Map
  gc = XCreateGC(dis,rw,0,0);

  for(i=0;i<127;i++){
    j=i*2+1;
    rgb[i+128].r=512+(int)j;
    rgb[i+128].g=512;
    rgb[i+128].b=512;
  }
  for(i=127;0<=i;i--){
    j=255-(i*2);
    rgb[i].r=512;
    rgb[i].g=512+(int)j;
    rgb[i].b=512;
  }

  rnm=0.0;
    for(k0=0;k0<istm0;k0=k0+1){
      for(y0=0;y0<YSIZE;y0++){
	for(x0=0;x0<XSIZE;x0++){
	  rr=dir[x0][y0][k0];
	  if(rr < 0.0){rr=-rr;}
	  if(rnm < rr){rnm=rr;}
	}
      }
    }
    //    printf(" max value %f \n", rnm);


    for(k0=0;k0<istm00;k0=k0+1){
      for(y0=0;y0<YSIZE;y0++){
	for(x0=0;x0<XSIZE;x0++){
	  dir_i=(int)((dir[x0][y0][k0]/rnm)*127+127);
	  DotRGB(dis,w,gc,x0*DOT_SIZE,WINDOW_HEIGHT-(y0+1)*DOT_SIZE,rgb[dir_i].r,rgb[dir_i].g,rgb[dir_i].b);
	}
      }
      usleep(8*10000);
    } 
    
  XFlush(dis);
  //  system("import -window Direction_Map rfs.bmp");
  
  printf("Hit Return Key\n");
  getchar();
  XFreeGC(dis,gc);
  XDestroyWindow(dis,w);
  XFlush(dis);
  XCloseDisplay(dis);
}

