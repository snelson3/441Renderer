#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
//#include "get_triangles.cxx"
#include <vtkPointData.h>
#include <math.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkDataSetWriter.h>

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


struct LightingParameters
{
    LightingParameters(void)
    {
         lightDir[0] = -0.6;
         lightDir[1] = 0;
         lightDir[2] = -0.8;
         Ka = 0.3;
         Kd = 0.7;
         Ks = 5.3;
         alpha = 7.5;
    };
  

    double lightDir[3]; // The direction of the light source
    double Ka;           // The coefficient for ambient lighting.
    double Kd;           // The coefficient for diffuse lighting.
    double Ks;           // The coefficient for specular lighting.
    double alpha;        // The exponent term for specular lighting.
};


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

struct Point
{
  double X;
  double Y;
  double Z;
  double color[3];
  double normal[3];
};


class Triangle
{
  public:
      double         X[3];
      double         Y[3];
      double         Z[3];
      double        colors[3][3];
      double        normals[3][3];

      void setPT1(double x, double y, double z){X[0]=x; Y[0]=y; Z[0] = z;};
      void setPT2(double x, double y, double z){X[1]=x; Y[1]=y; Z[1] = z;};
      void setPT3(double x, double y, double z){X[2]=x; Y[2]=y; Z[2] = z;};

      Point getPT(int t)
      {
        Point tr;
        tr.X = X[t];
        tr.Y = Y[t];
        tr.Z = Z[t];
        tr.color[0] = colors[t][0];
        tr.color[1] = colors[t][1];
        tr.color[2] = colors[t][2];
        tr.normal[0] = normals[t][0];//might be wrong, don't think so
        tr.normal[1] = normals[t][1];
        tr.normal[2] = normals[t][2];
        return tr;
      }

      double getX1(){return X[0];};
      double getX2(){return X[1];};
      double getX3(){return X[2];};

      double getY1(){return Y[0];};
      double getY2(){return Y[1];};
      double getY3(){return Y[2];};

      double getZ1(){return Z[0];};
      double getZ2(){return Z[1];};
      double getZ3(){return Z[2];};

      void setRed(double c,int i){colors[i][0] = c;};
      void setGreen(double c, int i){colors[i][1] = c;};
      void setBlue(double c, int i){colors[i][2] = c;};

      double getRed(int i){return colors[i][0];};
      double getGreen(int i){return colors[i][1];};
      double getBlue(int i){return colors[i][2];};

      void debugOut()
      {
            std::cerr<<"********************************"<<endl;
            std::cerr<<"V1("<<this->getX1()<<","<<this->getY1()<<","<<this->getZ1()<<")"<<endl;
            std::cerr<<"V2("<<this->getX2()<<","<<this->getY2()<<","<<this->getZ2()<<")"<<endl;
            std::cerr<<"V3("<<this->getX3()<<","<<this->getY3()<<","<<this->getZ3()<<")"<<endl;
      }

      void setV(Point *v1, Point *v2, Point *v3)
      {

        *v1 = getPT(0);
        *v2 = getPT(1);
        *v3 = getPT(2);
      }

      void sortVertices(Point *v1, Point *v2, Point *v3)
      {
        double vt1;
        double vt3;
        //could I use a bunch of for loops to reduce the logic here?

        vt1 = fmax(Y[0],fmax(Y[1],Y[2]));
        vt3 = fmin(Y[0],fmin(Y[1],Y[2]));
        if (vt1 == Y[0] && vt3 == Y[1])
          setV(v1,v3,v2);
        else if (vt1 == Y[0] && vt3 == Y[2])
          setV(v1,v2,v3);
        else if (vt1 == Y[1] && vt3 == Y[0])
          setV(v3,v1,v2);
        else if (vt1 == Y[1] && vt3 == Y[2])
          setV(v2,v1,v3);
        else if (vt1 == Y[2] && vt3 == Y[0])
          setV(v3,v2,v1);
        else if (vt1 == Y[2] && vt3 == Y[1])
          setV(v2,v3,v1);
        else
          std::cerr<<"SOMETHINGS BROkeN"<<endl;
      }
};

double dotProduct(double *A, double *B)
{
  return A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
}

std::vector<Triangle>
GetTriangles(void)
{
    vtkPolyDataReader *rdr = vtkPolyDataReader::New();
    rdr->SetFileName("proj1e_geometry.vtk");
    cerr << "Reading" << endl;
    rdr->Update();
    cerr << "Done reading" << endl;
    if (rdr->GetOutput()->GetNumberOfCells() == 0)
    {
        cerr << "Unable to open file!!" << endl;
        exit(EXIT_FAILURE);
    }
    vtkPolyData *pd = rdr->GetOutput();

    int numTris = pd->GetNumberOfCells();
    vtkPoints *pts = pd->GetPoints();
    vtkCellArray *cells = pd->GetPolys();
    vtkDoubleArray *var = (vtkDoubleArray *) pd->GetPointData()->GetArray("hardyglobal");
    double *color_ptr = var->GetPointer(0);
    //vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
    //float *color_ptr = var->GetPointer(0);
    vtkFloatArray *n = (vtkFloatArray *) pd->GetPointData()->GetNormals();
    float *normals = n->GetPointer(0);
    std::vector<Triangle> tris(numTris);
    vtkIdType npts;
    vtkIdType *ptIds;
    int idx;
    for (idx = 0, cells->InitTraversal() ; cells->GetNextCell(npts, ptIds) ; idx++)
    {
        if (npts != 3)
        {
            cerr << "Non-triangles!! ???" << endl;
            exit(EXIT_FAILURE);
        }
        double *pt = NULL;
        pt = pts->GetPoint(ptIds[0]);
        tris[idx].X[0] = (pt[0]+10)*50.0;
        tris[idx].Y[0] = (pt[1]+10)*50.0;
        tris[idx].Z[0] = (pt[2]-10)*0.05;
        tris[idx].normals[0][0] = normals[3*ptIds[0]+0];
        tris[idx].normals[0][1] = normals[3*ptIds[0]+1];
        tris[idx].normals[0][2] = normals[3*ptIds[0]+2];
        pt = pts->GetPoint(ptIds[1]);
        tris[idx].X[1] = (pt[0]+10)*50.0;
        tris[idx].Y[1] = (pt[1]+10)*50.0;
        tris[idx].Z[1] = (pt[2]-10)*0.05;
        tris[idx].normals[1][0] = normals[3*ptIds[1]+0];
        tris[idx].normals[1][1] = normals[3*ptIds[1]+1];
        tris[idx].normals[1][2] = normals[3*ptIds[1]+2];
        pt = pts->GetPoint(ptIds[2]);
        tris[idx].X[2] = (pt[0]+10)*50.0;
        tris[idx].Y[2] = (pt[1]+10)*50.0;
        tris[idx].Z[2] = (pt[2]-10)*0.05;
        tris[idx].normals[2][0] = normals[3*ptIds[2]+0];
        tris[idx].normals[2][1] = normals[3*ptIds[2]+1];
        tris[idx].normals[2][2] = normals[3*ptIds[2]+2];

        // 1->2 interpolate between light blue, dark blue
        // 2->2.5 interpolate between dark blue, cyan
        // 2.5->3 interpolate between cyan, green
        // 3->3.5 interpolate between green, yellow
        // 3.5->4 interpolate between yellow, orange
        // 4->5 interpolate between orange, brick
        // 5->6 interpolate between brick, salmon
        double mins[7] = { 1, 2, 2.5, 3, 3.5, 4, 5 };
        double maxs[7] = { 2, 2.5, 3, 3.5, 4, 5, 6 };
        unsigned char RGB[8][3] = { { 71, 71, 219 }, 
                                    { 0, 0, 91 },
                                    { 0, 255, 255 },
                                    { 0, 128, 0 },
                                    { 255, 255, 0 },
                                    { 255, 96, 0 },
                                    { 107, 0, 0 },
                                    { 224, 76, 76 } 
                                  };
        for (int j = 0 ; j < 3 ; j++)
        {
            float val = color_ptr[ptIds[j]];
            int r;
            for (r = 0 ; r < 7 ; r++)
            {
                if (mins[r] <= val && val < maxs[r])
                    break;
            }
            if (r == 7)
            {
                cerr << "Could not interpolate color for " << val << endl;
                exit(EXIT_FAILURE);
            }
            double proportion = (val-mins[r]) / (maxs[r]-mins[r]);
            tris[idx].colors[j][0] = (RGB[r][0]+proportion*(RGB[r+1][0]-RGB[r][0]))/255.0;
            tris[idx].colors[j][1] = (RGB[r][1]+proportion*(RGB[r+1][1]-RGB[r][1]))/255.0;
            tris[idx].colors[j][2] = (RGB[r][2]+proportion*(RGB[r+1][2]-RGB[r][2]))/255.0;
        }
    }

    return tris;
}


void debugOt(Point v1, Point v2, Point v3, std::vector<Triangle>::iterator i)
{
  cerr<<"***************"<<endl;
  cerr<<"v1("<<v1.X<<","<<v1.Y<<","<<v1.Z<<")"<<endl;
  cerr<<"v2("<<v2.X<<","<<v2.Y<<","<<v2.Z<<")"<<endl;
  cerr<<"v3("<<v3.X<<","<<v3.Y<<","<<v3.Z<<")"<<endl;
  cerr<<"**************"<<endl;
}


double findx (Point v1, Point v2, double E)
{
  return ((E-v1.Y)*(v2.X-v1.X))/(v2.Y-v1.Y)+v1.X;
}

double findc (double v1c, double v2c, double v1x, double v2x, double F)
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
        //View direction is going to be [0,0,-1];
        double V[3];
        V[0] = 0;
        V[1] = 0;
        V[2] = -1;
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
        return shading;
      }

      void setColor(double x, double y, Point v1, Point v2, Point v3)
      {
        if (x == this->width)
            return;
        double v5x = findx(v1,v2,y);
        double v6x = findx(v2,v3,y);

        double v5z = findc(v1.Z,v2.Z,v1.Y,v2.Y,y);
        double v6z = findc(v2.Z,v3.Z,v2.Y,v3.Y,y);
        double v4z = findc(v5z,v6z,v5x,v6x,x);

         if ((v4z <= getZ(x,y) ))
            return;

        setZ(x,y,v4z);

        double v5r = findc(v1.color[0],v2.color[0],v1.Y,v2.Y,y);
        double v6r = findc(v2.color[0],v3.color[0],v2.Y,v3.Y,y);
        double v4r = findc(v5r,v6r,v5x,v6x,x);

        double v5u = findc(v1.color[2],v2.color[2],v1.Y,v2.Y,y);
        double v6u = findc(v2.color[2],v3.color[2],v2.Y,v3.Y,y);
        double v4u = findc(v5u,v6u,v5x,v6x,x);

        double v5g = findc(v1.color[1],v2.color[1],v1.Y,v2.Y,y);
        double v6g = findc(v2.color[1],v3.color[1],v2.Y,v3.Y,y);
        double v4g = findc(v5g,v6g,v5x,v6x,x);

        double v5n[3];
        double v6n[3];
        double v4n[3];

        v5n[0] = findc(v1.normal[0],v2.normal[0],v1.Y,v2.Y,y);
        v5n[1] = findc(v1.normal[1],v2.normal[1],v1.Y,v2.Y,y);
        v5n[2] = findc(v1.normal[2],v2.normal[2],v1.Y,v2.Y,y);

        v6n[0] = findc(v2.normal[0],v3.normal[0],v2.Y,v3.Y,y);
        v6n[1] = findc(v2.normal[1],v3.normal[1],v2.Y,v3.Y,y);
        v6n[2] = findc(v2.normal[2],v3.normal[2],v2.Y,v3.Y,y);

        v4n[0] = findc(v5n[0],v6n[0],v5x,v6x,x);
        v4n[1] = findc(v5n[1],v6n[1],v5x,v6x,x);
        v4n[2] = findc(v5n[2],v6n[2],v5x,v6x,x);

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
          if (x >= 0 && x< screen.width)
            screen.setColor(x,y,v1,v2,v3);
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


void depositTriangle(std::vector<Triangle>::iterator i, Screen screen)
{
  Point v1,v2,v3,v4;
  double slope, intercept;

  /*Step 1*/

  //3 Cases, when 1 and 2 are flat bottom and when 1 and 3 are flat bottom and when 2 and 3 are flat bottom
  if (i->getY1() == i->getY2())
  {
    if (i->getX1() < i->getX2())
      i->setV(&v1,&v3,&v2);
    else
      i->setV(&v3,&v1,&v2);

    if (v2.Y < v1.Y)
      writeFlatTop(v1,v2,v3,screen);
    else
      writeFlatBottom(v1,v2,v3,screen);
  }
  else if (i->getY2() == i->getY3())
  {
    if (i->getX2() < i->getX3())
      i->setV(&v2,&v1,&v3);
    else
      i->setV(&v2,&v3,&v1);
    
    if (v2.Y < v1.Y)
      writeFlatTop(v1,v2,v3,screen);
    else
      writeFlatBottom(v1,v2,v3,screen);
  }
  else if (i->getY1() == i->getY3())
  {
    if (i->getX1() < i->getX3())
      i->setV(&v1,&v2,&v3);
    else
      i->setV(&v3,&v2,&v1);

    if (v2.Y < v1.Y)
      writeFlatTop(v1,v2,v3,screen);
    else
      writeFlatBottom(v1,v2,v3,screen);
  }
  else
  {
    //must split into a flat bottom and a flat top triangle

    /*Step 2*/


    i->sortVertices(&v1,&v2,&v3);

    /*Step 3*/

    if (v1.X == v3.X) slope = 0; //right triangle
    else slope = (v1.Y - v3.Y)/(v1.X - v3.X);
    intercept = v3.Y - v3.X*slope;

    /*Step 4*/

    v4.Y = v2.Y;

    if (slope == 0) v4.X = v3.X;
    else v4.X = (v4.Y - intercept) / slope;

    /*Step 4.5*/
    
    v4.Z = findc(v1.Z,v3.Z,v1.Y,v3.Y,v4.Y);

    v4.color[0] = findc(v1.color[0],v3.color[0],v1.Y,v3.Y,v4.Y);
    v4.color[1] = findc(v1.color[1],v3.color[1],v1.Y,v3.Y,v4.Y);
    v4.color[2] = findc(v1.color[2],v3.color[2],v1.Y,v3.Y,v4.Y);

    v4.normal[0] = findc(v1.normal[0],v3.normal[0],v1.Y,v3.Y,v4.Y);
    v4.normal[1] = findc(v1.normal[1],v3.normal[1],v1.Y,v3.Y,v4.Y);
    v4.normal[2] = findc(v1.normal[2],v3.normal[2],v1.Y,v3.Y,v4.Y);          
    
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

int main()
{
   vtkImageData *image = NewImage(1000, 1000);
   unsigned char *buffer =
     (unsigned char *) image->GetScalarPointer(0,0,0);

   int npixels = 1000*1000;
    double *zbuffer = new double[npixels];
   // for (int i = 0; i < npixels*3 ; i++)
   //     buffer[i] = 0; //You can use loop unrolling so this is not 2 for loops
   // for (int i = 0; i < npixels; i++)
   //    zbuffer[i] = -1; //valid depth values go from -1 (back) to 0 (front)

    for (int i = 0; i <npixels*3; i+=3)
    {
        buffer[i] = 0;
        buffer[i+1] = 0;
        buffer[i+2] = 0;
        zbuffer[i] = -1;
    }

   std::vector<Triangle> triangles = GetTriangles();

   Screen screen;
   screen.buffer = buffer;
   screen.zbuffer = zbuffer;
   screen.width = 1000;
   screen.height = 1000;

   for(std::vector<Triangle>::iterator i = triangles.begin(); i != triangles.end(); ++i)
     depositTriangle(i,screen);
   
   WriteImage(image, "allTriangles");
}
