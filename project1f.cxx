#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <math.h>
#include "reader1F.cxx"
#include "shading.cxx"
#include "matrix.cxx"

using std::cerr;
using std::endl;

bool k = false;

double ceil441(double f)
{
    return ceil(f-0.00001);
}

double floor441(double f)
{
    return floor(f+0.00001);
}

vtkImageData *
NewImage(int width, int height)
{
    vtkImageData *img = vtkImageData::New();
    img->SetDimensions(width, height, 1);
    img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    return img;
}

void
WriteImage(vtkImageData *img, const char *filename)
{
   std::string full_filename = filename;
   full_filename += ".png";
   vtkPNGWriter *writer = vtkPNGWriter::New();
   writer->SetInputData(img);
   writer->SetFileName(full_filename.c_str());
   writer->Write();
   writer->Delete();
}

double cotan(double deg)
{
  return 1/tan(deg);
}

double dotProduct(double *A, double *B)
{
  return A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
}

void crossProduct(double *A, double *B, double *AB)
{
  AB[0] = A[1]*B[2]-A[2]*B[1];
  AB[1] = B[0]*A[2]-A[0]*B[2];
  AB[2] = A[0]*B[1]-A[1]*B[0];
}

void debugOt(Point v1, Point v2, Point v3, std::vector<Triangle>::iterator i)
{
  cerr<<"***************"<<endl;
  cerr<<"v1("<<v1.X<<","<<v1.Y<<","<<v1.Z<<")"<<endl;
  cerr<<"v2("<<v2.X<<","<<v2.Y<<","<<v2.Z<<")"<<endl;
  cerr<<"v3("<<v3.X<<","<<v3.Y<<","<<v3.Z<<")"<<endl;
  cerr<<"**************"<<endl;
}

void normalize(double *A)
{
  double norm = sqrt(A[0]*A[0]+A[1]*A[1]+A[2]*A[2]);
  A[0] = A[0]/norm;
  A[1] = A[1]/norm;
  A[2] = A[2]/norm;
}

double findx (Point v1, Point v2, double E)
{
  return ((E-v1.Y)*(v2.X-v1.X))/(v2.Y-v1.Y)+v1.X;
}

double interpolate (double v1c, double v2c, double v1x, double v2x, double F)
{
  return v1c+((F-v1x)/(v2x-v1x))*(v2c-v1c);
}

class Screen
{
  public:
      unsigned char   *buffer;
      double *zbuffer;
      int width, height;
      LightingParameters lp;

      void setRed(double r, int x, int y)
      {
        int loc = y * (width)+x;
        buffer[loc*3] = (unsigned char)ceil441(fmin(1,r)*255.0);
      }

      void setGreen(double g, int x, int y)
      {
        int loc = y * (width) + x;
        buffer[(loc*3)+1] = (unsigned char)ceil441(fmin(1,g)*255.0);
      }
      void setBlue(double u, int x, int y){
        int loc = y * (width) + x;
        buffer[(loc*3)+2] = (unsigned char)ceil441(fmin(1,u)*255.0);
      }

      double getZ(int x, int y){
        int loc = y * (width) + x;
        return zbuffer[loc];
      }

      void setZ(int x, int y, double z){
        int loc = y *(width) + x;
        zbuffer[loc] = z;
      }

      double calculateDiffuse(double *normal)
      {
        double diff = fabs(dotProduct(lp.lightDir,normal));
        return diff;
      }

      double calculateSpecular(double *normal)
      {
        //View direction is going to be [0,0,-1];//jk not really
        double V[3];
        V[0] = 0;
        V[1] = 0;
        V[2] = 1;
        double R[3];
        R[0] = 2 * dotProduct(lp.lightDir,normal)*normal[0]-lp.lightDir[0];
        R[1] = 2 * dotProduct(lp.lightDir,normal)*normal[1]-lp.lightDir[1];
        R[2] = 2 * dotProduct(lp.lightDir,normal)*normal[2]-lp.lightDir[2];
        return dotProduct(R,V);
      }

      double calculateShading(double x, double y, double *normal)
      {
        double diffuse = calculateDiffuse(normal);
        double specular = fmax(0,pow(calculateSpecular(normal),lp.alpha));
        double shading = lp.Ka + lp.Kd*diffuse + lp.Ks*specular;
        return 1;//shading;
      }

      void setColor(double x, double y, Point v1, Point v2, Point v3)
      {
        if (x == this->width)
            return;
        double v5x = findx(v1,v2,y);
        double v6x = findx(v2,v3,y);

        double v5z = interpolate(v1.Z,v2.Z,v1.Y,v2.Y,y);
        double v6z = interpolate(v2.Z,v3.Z,v2.Y,v3.Y,y);
        double v4z = interpolate(v5z,v6z,v5x,v6x,x);

        if ((v4z <= getZ(x,y) ))//I think I'm supposed to reverse this (I already did, it was originally <=)
             return;

        setZ(x,y,v4z);

        double v5r = interpolate(v1.color[0],v2.color[0],v1.Y,v2.Y,y);
        double v6r = interpolate(v2.color[0],v3.color[0],v2.Y,v3.Y,y);
        double v4r = interpolate(v5r,v6r,v5x,v6x,x);

        double v5u = interpolate(v1.color[2],v2.color[2],v1.Y,v2.Y,y);
        double v6u = interpolate(v2.color[2],v3.color[2],v2.Y,v3.Y,y);
        double v4u = interpolate(v5u,v6u,v5x,v6x,x);

        double v5g = interpolate(v1.color[1],v2.color[1],v1.Y,v2.Y,y);
        double v6g = interpolate(v2.color[1],v3.color[1],v2.Y,v3.Y,y);
        double v4g = interpolate(v5g,v6g,v5x,v6x,x);

        double v5n[3];
        double v6n[3];
        double v4n[3];

        v5n[0] = interpolate(v1.normal[0],v2.normal[0],v1.Y,v2.Y,y);
        v5n[1] = interpolate(v1.normal[1],v2.normal[1],v1.Y,v2.Y,y);
        v5n[2] = interpolate(v1.normal[2],v2.normal[2],v1.Y,v2.Y,y);

        v6n[0] = interpolate(v2.normal[0],v3.normal[0],v2.Y,v3.Y,y);
        v6n[1] = interpolate(v2.normal[1],v3.normal[1],v2.Y,v3.Y,y);
        v6n[2] = interpolate(v2.normal[2],v3.normal[2],v2.Y,v3.Y,y);

        v4n[0] = interpolate(v5n[0],v6n[0],v5x,v6x,x);
        v4n[1] = interpolate(v5n[1],v6n[1],v5x,v6x,x);
        v4n[2] = interpolate(v5n[2],v6n[2],v5x,v6x,x);

        this->setRed(v4r*this->calculateShading(x,y,v4n), x, y);
        this->setBlue(v4u*this->calculateShading(x,y,v4n), x, y);
        this->setGreen(v4g*this->calculateShading(x,y,v4n), x, y);
        }
};

void writeFlatBottom(Point v1, Point v2, Point v3, Screen screen)
{
    double leftslope, rightslope, li, ri, leftend, rightend;

    /********FIND THE LEFT AND RIGHT SLOPE AND INTERCEPT***********/
    if (v1.X==v2.X)
      leftslope = 0;//right triangle
    else
      leftslope = (v2.Y-v1.Y)/(v2.X-v1.X);
    li = v1.Y-leftslope*v1.X;

    if (v2.X == v3.X)
      rightslope = 0; //right triangle
    else
      rightslope = (v2.Y - v3.Y)/(v2.X - v3.X);
    ri = v3.Y - v3.X*rightslope;

    /****Loop through every hor line of the triangle***/

    for (int y = ceil441(v1.Y); y <= floor441(v2.Y); y++)
    {
      if (y >= 0 && y < screen.height)
      {
        //now I need to find the left endpoint and the right endpoint
        if (leftslope == 0) leftend = v1.X;
        else leftend = (y - li) / leftslope;
        if (rightslope == 0) rightend = v3.X;
        else rightend = (y - ri) / rightslope;

        if (ceil441(leftend) == floor(rightend) && (leftend >= 0 && leftend < screen.width))
            screen.setColor(ceil(leftend),y,v1,v2,v3);

        for (int x = ceil441(leftend); x <= floor441(rightend); x++)
        {
         // cerr<<"WEE";
          if (x >= 0 && x< screen.width)
            screen.setColor(x,y,v1,v2,v3);
        }
      }
    }
}

void writeFlatTop(Point v1, Point v2, Point v3, Screen screen)
{
    double leftslope, rightslope, li, ri,leftend,rightend;

    if (v1.X==v2.X)
      leftslope = 0;//right triangle
    else
      leftslope = (v2.Y-v1.Y)/(v2.X-v1.X);
    li = v1.Y-leftslope*v1.X;

    if (v2.X == v3.X)
      rightslope = 0; //right triangle
    else
      rightslope = (v2.Y - v3.Y)/(v2.X - v3.X);
    ri = v3.Y - v3.X * rightslope;

    /****Loop through every hor line of the Triangle***/

    for (int y = floor441(v1.Y); y >= ceil441(v2.Y); y--)
    {
      if (y >= 0 && y < screen.height)
      {
      //now I need to find the left endpoint and the right endpoint
        if (leftslope == 0) leftend = v1.X;
        else leftend = (y - li) / leftslope;
        if (rightslope == 0) rightend = v3.X;
        else rightend = (y - ri) / rightslope;

        if ( ceil441(leftend) == floor441(rightend) && ( (leftend >= 0) && (leftend< screen.width) ) )
          screen.setColor(ceil441(leftend),y,v1,v2,v3);

        for (int x = ceil441(leftend); x <= floor441(rightend); x++)
          if (x >= 0 && x< screen.width)
            screen.setColor(x,y,v1,v2,v3);
      }
    }
}


void depositTriangle(Triangle i, Screen screen)
{
  Point v1,v2,v3,v4;
  double slope, intercept;

  /*Step 1*/

  //3 Cases, when 1 and 2 are flat bottom and when 1 and 3 are flat bottom and when 2 and 3 are flat bottom
  if (i.getY1() == i.getY2())
  {
    if (i.getX1() < i.getX2())
      i.setV(&v1,&v3,&v2);
    else
      i.setV(&v3,&v1,&v2);

    if (v2.Y < v1.Y)
      writeFlatTop(v1,v2,v3,screen);
    else
      writeFlatBottom(v1,v2,v3,screen);
  }
  else if (i.getY2() == i.getY3())
  {
    if (i.getX2() < i.getX3())
      i.setV(&v2,&v1,&v3);
    else
      i.setV(&v2,&v3,&v1);
    
    if (v2.Y < v1.Y)
      writeFlatTop(v1,v2,v3,screen);
    else
      writeFlatBottom(v1,v2,v3,screen);
  }
  else if (i.getY1() == i.getY3())
  {
    if (i.getX1() < i.getX3())
      i.setV(&v1,&v2,&v3);
    else
      i.setV(&v3,&v2,&v1);

    if (v2.Y < v1.Y)
      writeFlatTop(v1,v2,v3,screen);
    else
      writeFlatBottom(v1,v2,v3,screen);
  }
  else
  {
    //must split into a flat bottom and a flat top triangle

    /*Step 2*/


    i.sortVertices(&v1,&v2,&v3);

    /*Step 3*/

    if (v1.X == v3.X) slope = 0; //right triangle
    else slope = (v1.Y - v3.Y)/(v1.X - v3.X);
    intercept = v3.Y - v3.X*slope;

    /*Step 4*/

    v4.Y = v2.Y;

    if (slope == 0) v4.X = v3.X;
    else v4.X = (v4.Y - intercept) / slope;

    /*Step 4.5*/
    
    v4.Z = interpolate(v1.Z,v3.Z,v1.Y,v3.Y,v4.Y);

    v4.color[0] = interpolate(v1.color[0],v3.color[0],v1.Y,v3.Y,v4.Y);
    v4.color[1] = interpolate(v1.color[1],v3.color[1],v1.Y,v3.Y,v4.Y);
    v4.color[2] = interpolate(v1.color[2],v3.color[2],v1.Y,v3.Y,v4.Y);

    v4.normal[0] = interpolate(v1.normal[0],v3.normal[0],v1.Y,v3.Y,v4.Y);
    v4.normal[1] = interpolate(v1.normal[1],v3.normal[1],v1.Y,v3.Y,v4.Y);
    v4.normal[2] = interpolate(v1.normal[2],v3.normal[2],v1.Y,v3.Y,v4.Y);          
    
    /*Step 5*/

    if (v2.X < v4.X)
      writeFlatTop(v2,v3,v4,screen);
    else
      writeFlatTop(v4,v3,v2,screen);

    /*Step 6*/

    if (v2.X<v4.X)
      writeFlatBottom(v2,v1,v4,screen);
    else 
      writeFlatBottom(v4,v1,v2,screen);
  }
}



class Camera
{
  public:
    double          near, far;
    double          angle;
    double          position[3];
    double          focus[3];
    double          up[3];
    Matrix          ct;
    Matrix          vt;
    Matrix          dt;
    Matrix          tt;


    void          ViewTransform(void) 
    {
      vt.A[0][0] = cotan(angle/2);
      vt.A[0][1] = 0;
      vt.A[0][2] = 0;
      vt.A[0][3] = 0;

      vt.A[1][0] = 0;
      vt.A[1][1] = cotan(angle/2);
      vt.A[1][2] = 0;
      vt.A[1][3] = 0;

      vt.A[2][0] = 0;
      vt.A[2][1] = 0;
      vt.A[2][2] = (far+near)/(far-near);
      vt.A[2][3] = -1;

      vt.A[3][0] = 0;
      vt.A[3][1] = 0;
      vt.A[3][2] = (2*far*near)/(far-near);
      vt.A[3][3] = 0;
    };


    void          CameraTransform(double *v1, double *v2, double *v3, double *o) 
    {
      ct.A[0][0] = v1[0];
      ct.A[0][1] = v2[0];
      ct.A[0][2] = v3[0];
      ct.A[0][3] = 0;

      ct.A[1][0] = v1[1];
      ct.A[1][1] = v2[1];
      ct.A[1][2] = v3[1];
      ct.A[1][3] = 0;

      ct.A[2][0] = v1[2];
      ct.A[2][1] = v2[2];
      ct.A[2][2] = v3[2];
      ct.A[2][3] = 0;

      double t[3];
      for (int i = 0; i < 3; i++) t[i] = 0-o[i];

      ct.A[3][0] = dotProduct(v1,t);
      ct.A[3][1] = dotProduct(v2,t);
      ct.A[3][2] = dotProduct(v3,t);
      ct.A[3][3] = 1;

    };



    void          DeviceTransform(double *A) 
    {

      double x = (1000*(A[0]+1))/2;
      double y = (1000*(A[1]+1))/2;
      double z = A[2];
      dt.A[0][0] = 1000/2;
      dt.A[0][1] = 0;
      dt.A[0][2] = 0;
      dt.A[0][3] = 0;

      dt.A[1][0] = 0;
      dt.A[1][1] = 1000/2;
      dt.A[1][2] = 0;
      dt.A[1][3] = 0;

      dt.A[2][0] = 0;
      dt.A[2][1] = 0;
      dt.A[2][2] = 1;
      dt.A[2][3] = 0;

      dt.A[3][0] = 1000/2;
      dt.A[3][1] = 1000/2;
      dt.A[3][2] = 0;
      dt.A[3][3] = 1;
    };
};

//calculate the camera frame somehow

double SineParameterize(int curFrame, int nFrames, int ramp)
{
    int nNonRamp = nFrames-2*ramp;
    double height = 1./(nNonRamp + 4*ramp/M_PI);
    if (curFrame < ramp)
    {
        double factor = 2*height*ramp/M_PI;
        double eval = cos(M_PI/2*((double)curFrame)/ramp);
        return (1.-eval)*factor;
    }
    else if (curFrame > nFrames-ramp)
    {
        int amount_left = nFrames-curFrame;
        double factor = 2*height*ramp/M_PI;
        double eval =cos(M_PI/2*((double)amount_left/ramp));
        return 1. - (1-eval)*factor;
    }
    double amount_in_quad = ((double)curFrame-ramp);
    double quad_part = amount_in_quad*height;
    double curve_part = height*(2*ramp)/M_PI;
    return quad_part+curve_part;
}

Camera
GetCamera(int frame, int nframes)
{
    double t = SineParameterize(frame, nframes, nframes/10);
    Camera c;
    c.near = 5;
    c.far = 200;
    c.angle = M_PI/6;
    c.position[0] = 40*sin(2*M_PI*t);
    c.position[1] = 40*cos(2*M_PI*t);
    c.position[2] = 40;
    c.focus[0] = 0;
    c.focus[1] = 0;
    c.focus[2] = 0;
    c.up[0] = 0;
    c.up[1] = 1;
    c.up[2] = 0;
    // cerr<<"c.near "<<c.near<<endl;
    // cerr<<"c.far "<<c.far<<endl;
    // cerr<<"c.angle "<<c.angle<<endl;

    double          v1[3];
    double          v2[3];
    double          v3[3];
        //Up x (O-focus)      V1
        //(O-focus) x v1      V2
        //(O-focus)           V3

    for (int i = 0; i < 3; i++) v3[i] = c.position[i] - c.focus[i];
    crossProduct(c.up,v3,v1);
    normalize(v1);
    crossProduct(v3,v1,v2);
    normalize(v2);
    normalize(v3);

    // cerr<<"v1 " << v1[0] << " "<<v1[1]<<" "<<v1[2]<<endl;
    // cerr<<"v2 " << v2[0] << " "<<v2[1]<<" "<<v2[2]<<endl;
    // cerr<<"v3 " << v3[0] << " "<<v3[1]<<" "<<v3[2]<<endl;
    // cerr<<"o " << c.position[0] << " "<<c.position[1]<<" "<<c.position[2]<<endl;



    //dt ct vt
    //dt vt ct
    c.CameraTransform(v1,v2,v3,c.position);
    c.ViewTransform();
    c.DeviceTransform(c.position);
    c.tt = c.tt.ComposeMatrices(c.ct,c.vt);
    // cerr<<"*********Combined ct with vt**************"<<endl;
    // c.tt.Print(cerr);
          c.tt = c.tt.ComposeMatrices(c.tt,c.dt);
    // cerr<<"*********Combined tt with dt***************"<<endl;
    // c.tt.Print(cerr);
    //c.vt.Print(cerr);
    // camtr.Print(cerr);
    return c;
}






Triangle toDeviceSpace(std::vector<Triangle>::iterator i, Camera c)
  {
    Triangle t;

    double oPT1[4];
    oPT1[0]=i->getX1(); oPT1[1]=i->getY1(); oPT1[2]=i->getZ1(); oPT1[3]=1;
    double oPT2[4];
    oPT2[0]=i->getX2(); oPT2[1]=i->getY2(); oPT2[2]=i->getZ2(); oPT2[3]=1;
    double oPT3[4];
    oPT3[0]=i->getX3(); oPT3[1]=i->getY3(); oPT3[2]=i->getZ3(); oPT3[3]=1;

    double nPT1[4];
    double nPT2[4];
    double nPT3[4];
    for (int i = 0; i<4; i++)
    {
      nPT1[i] = 0; nPT2[i] = 0; nPT3[i] = 0;
    }

    //c.tt.Print(cerr);


    c.tt.TransformPoint(oPT1,nPT1);
    c.tt.TransformPoint(oPT2,nPT2);
    c.tt.TransformPoint(oPT3,nPT3);

    double d = 1;
    double m = 0;

    // normalize(nPT1);
    // normalize(nPT2);
    // normalize(nPT3);

    t.setPT1(nPT1[0]/nPT1[3],nPT1[1]/nPT1[3],nPT1[2]/nPT1[3]);
    t.setPT2(nPT2[0]/nPT2[3],nPT2[1]/nPT1[3],nPT2[2]/nPT1[3]);
    t.setPT3(nPT3[0]/nPT3[3],nPT3[1]/nPT1[3],nPT3[2]/nPT1[3]);
    i->setColors(t.colors);
    i->setNormals(t.normals);


    // double T[4];
    // T[0] = 1.11111;
    // T[1] = 7.46665;
    // T[2] = -8.9899;
    // T[3] = 1;

    // double nT[4];

    // c.tt.TransformPoint(T,nT);

    // cerr<<"Transformed " <<T[0] <<" "<<T[1] <<" "<<T[2]<<" "<<T[3]<<" to "<<nT[0]/nT[3]<<" "<<nT[1]/nT[3]<<" "<<nT[2]/nT[3]<<" "
    //     <<nT[3]<<endl;

    return t;
  }

void makeImage(int curr, int total, std::vector<Triangle> t)
{
  //set up the screen
vtkImageData *image = NewImage(1000, 1000);
   unsigned char *buffer =
     (unsigned char *) image->GetScalarPointer(0,0,0);

   int npixels = 1000*1000;
    double *zbuffer = new double[npixels];
    for (int i = 0; i < npixels*3 ; i++)
        buffer[i] = 0;
    for (int i = 0; i < npixels; i++)
        zbuffer[i] = -1; //I was told depth is opposite in this part

   Screen screen;
   screen.buffer = buffer;
   screen.zbuffer = zbuffer;
   screen.width = 1000;
   screen.height = 1000;

  //Now to get the image, put each triangle in device space, and render it on the screen (before writing the png)

  Camera c = GetCamera(curr,total);
  //I probably need to do something to get the right Matrix here to send into toDeviceSpace
  int e = 0;
  for(std::vector<Triangle>::iterator i = t.begin(); i != t.end(); ++i)
  {
    // if (e != 0)
    //   continue;
    // // e++;
    Triangle t = toDeviceSpace(i,c);
    depositTriangle(t,screen);
    if (t.getX1() > e)
      e = t.getX1();

     // cerr<<"************ t points **************"<<endl;
     // cerr<<t.getX1()<<" "<<t.getY1()<<" "<<t.getZ1()<<endl;
     // cerr<<t.getX2()<<" "<<t.getY2()<<" "<<t.getZ2()<<endl;
     // cerr<<t.getX3()<<" "<<t.getY3()<<" "<<t.getZ3()<<endl;
    // cerr<<"************ t colors **************"<<endl;
    // cerr<<t.colors[0][0]<<" "<<t.colors[0][1]<<" "<<t.colors[0][2]<<endl;
    // cerr<<t.colors[1][0]<<" "<<t.colors[1][1]<<" "<<t.colors[1][2]<<endl;
    // cerr<<t.colors[2][0]<<" "<<t.colors[2][1]<<" "<<t.colors[2][2]<<endl;
    // cerr<<"************ t normals **************"<<endl;
    // cerr<<t.normals[0][0]<<" "<<t.normals[0][1]<<" "<<t.normals[0][2]<<endl;
    // cerr<<t.normals[1][0]<<" "<<t.normals[1][1]<<" "<<t.normals[1][2]<<endl;
    // cerr<<t.normals[2][0]<<" "<<t.normals[2][1]<<" "<<t.normals[2][2]<<endl;
    // cerr<<"*************************************"<<endl;
    // cerr<<"************ i points **************"<<endl;
    // cerr<<i->getX1()<<" "<<i->getY1()<<" "<<i->getZ1()<<endl;
    // cerr<<i->getX2()<<" "<<i->getY2()<<" "<<i->getZ2()<<endl;
    // cerr<<i->getX3()<<" "<<i->getY3()<<" "<<i->getZ3()<<endl;
    // cerr<<"************ i colors **************"<<endl;
    // cerr<<i->colors[0][0]<<" "<<i->colors[0][1]<<" "<<i->colors[0][2]<<endl;
    // cerr<<i->colors[1][0]<<" "<<i->colors[1][1]<<" "<<i->colors[1][2]<<endl;
    // cerr<<i->colors[2][0]<<" "<<i->colors[2][1]<<" "<<i->colors[2][2]<<endl;
    // cerr<<"************ i normals **************"<<endl;
    // cerr<<i->normals[0][0]<<" "<<i->normals[0][1]<<" "<<i->normals[0][2]<<endl;
    // cerr<<i->normals[1][0]<<" "<<i->normals[1][1]<<" "<<i->normals[1][2]<<endl;
    // cerr<<i->normals[2][0]<<" "<<i->normals[2][1]<<" "<<i->normals[2][2]<<endl;
    // cerr<<"*************************************"<<endl;

  }
  cerr<<"E "<<e<<endl;
  //This breakfs with big amounts of frames but thats ok for now/this project
  char framenum[6];
  if (curr < 10) sprintf(framenum, "fr000%d",curr);
  else if (curr < 100) sprintf(framenum, "fr00%d",curr);
  else if (curr < 1000) sprintf(framenum, "fr0%d",curr);
  else  sprintf(framenum, "fr%d",curr);
  WriteImage(image, framenum);
}

int main()
{
  std::vector<Triangle> triangles = GetTriangles();
  makeImage(0,1000,triangles);
  
  // Camera c;
  // c.near = 5;
  // c.far = 10;
  // c.angle = 90;

  // c.ViewTransform();
  // c.vt.Print(cerr);

//  cerr<<cotan(90/2);

   makeImage(250,1000,triangles);
   makeImage(500,1000,triangles);
   makeImage(750,1000,triangles);


   // vtkImageData *image = NewImage(1000, 1000);
   // unsigned char *buffer =
   //   (unsigned char *) image->GetScalarPointer(0,0,0);

   // int npixels = 1000*1000;
   //  double *zbuffer = new double[npixels];
   //  for (int i = 0; i < npixels*3 ; i++)
   //      buffer[i] = 0;
   //  for (int i = 0; i < npixels; i++)
   //      zbuffer[i] = -1; //valid depth values go from -1 (back) to 0 (front)

   // Screen screen;
   // screen.buffer = buffer;
   // screen.zbuffer = zbuffer;
   // screen.width = 1000;
   // screen.height = 1000;

   // for(std::vector<Triangle>::iterator i = triangles.begin(); i != triangles.end(); ++i)
   //   depositTriangle(i,screen);
   
   // WriteImage(image, "allTriangles");
}
