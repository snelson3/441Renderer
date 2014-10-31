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
#include <vtkCellArray.h>


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

struct Point
{
  double X;
  double Y;
  double Z;
  bool isv4;
  double color[3];
};


class Triangle
{
  public:
      double         X[3];
      double         Y[3];
      double         Z[3];
      double        colors[3][3];

      // Point PT1;
      // Point PT2;
      // Point PT3;

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

      //This is not the correct places to put the colors I think
        //I'm not even sure why the pixels are a 2d array now.
      void setRed(double c,int i){colors[i][0] = c;};
      void setGreen(double c, int i){colors[i][1] = c;};
      void setBlue(double c, int i){colors[i][2] = c;};

      //these need to be changed too.
      double getRed(int i){return colors[i][0];};
      double getGreen(int i){return colors[i][1];};
      double getBlue(int i){return colors[i][2];};

      void debugOut()
      {
            std::cerr<<"********************************"<<endl;
            std::cerr<<"V1("<<this->getX1()<<","<<this->getY1()<<","<<this->getZ1()<<")"<<endl;
            std::cerr<<"V2("<<this->getX2()<<","<<this->getY2()<<","<<this->getZ2()<<")"<<endl;
            std::cerr<<"V3("<<this->getX3()<<","<<this->getY3()<<","<<this->getZ3()<<")"<<endl;
            //std::cerr<<"C("<<this->getRed()<<","<<this->getBlue()<<","<<this->getGreen()<<")"<<endl;
      }

      void setV(Point *v1, Point *v2, Point *v3)
      {

        *v1 = getPT(0);
        *v2 = getPT(1);
        *v3 = getPT(2);
        // PT1 = *v1;
        // PT2 = *v2;
        // PT3 = *v3;//I think its because
      }

      void sortVertices(Point *v1, Point *v2, Point *v3)
      {
        double vt1;
        double vt3;
        //could I use a bunch of for loops to reduce the logic here?

        vt1 = fmax(Y[0],fmax(Y[1],Y[2]));
        vt3 = fmin(Y[0],fmin(Y[1],Y[2]));
        if (vt1 == Y[0] && vt3 == Y[1])
        {
          setV(v1,v3,v2);

        }
        else if (vt1 == Y[0] && vt3 == Y[2])
        {
          setV(v1,v2,v3);
        }
        else if (vt1 == Y[1] && vt3 == Y[0])
        {
          setV(v3,v1,v2);
        }
        else if (vt1 == Y[1] && vt3 == Y[2])
        {
          setV(v2,v1,v3);
        }
        else if (vt1 == Y[2] && vt3 == Y[0])
        {
          setV(v3,v2,v1);
        }
        else if (vt1 == Y[2] && vt3 == Y[1])
        {
          setV(v2,v3,v1);
        }
        else
          std::cerr<<"SOMETHINGS BROkeN"<<endl;
      }

  // would some methods for transforming the triangle in place be helpful?
};

void debugOt(Point v1, Point v2, Point v3, std::vector<Triangle>::iterator i)
{
  cerr<<"***************"<<endl;
  cerr<<"v1("<<v1.X<<","<<v1.Y<<","<<v1.Z<<")"<<endl;
  cerr<<"v2("<<v2.X<<","<<v2.Y<<","<<v2.Z<<")"<<endl;
  cerr<<"v3("<<v3.X<<","<<v3.Y<<","<<v3.Z<<")"<<endl;
  cerr<<"**************"<<endl;
  // cerr<<"TRIANGLE OF HELL";
  // cerr<<endl<<"v1("<<i->PT1.X<<","<<i->PT1.Y<<","<<i->PT1.Z<<")"<<endl;
  // cerr<<"v2("<<i->PT2.X<<","<<i->PT2.Y<<","<<i->PT2.Z<<")"<<endl;
  // cerr<<"v3("<<i->PT3.X<<","<<i->PT3.Y<<","<<i->PT3.Z<<")"<<endl;
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

      //assume 9 pixels, so setr[3][3] would be spot 24

      void setRed(double r, int x, int y){
        //std::cerr<<"  SETRED "<<x<<endl;
        int loc = y * (width)+x;
        buffer[loc*3] = (unsigned char)ceil441(r*255.0);
        //buffer[(y*width)+x] = r;
      }
      void setGreen(double u, int x, int y){
        int loc = y * (width) + x;
        buffer[(loc*3)+1] = (unsigned char)ceil441(u*255.0);
        //buffer[(y*width)+x+1] = u;
      }
      void setBlue(double g, int x, int y){
        int loc = y * (width) + x;
        //buffer[(y*width)+x+2]=g;
        buffer[(loc*3)+2] = (unsigned char)ceil441(g*255.0);
        //cerr<<x<<y<<"\n";
      }

      double getZ(int x, int y){
        int loc = y * (width) + x;
        return zbuffer[loc];
      }

      void setZ(int x, int y, double z){
        int loc = y *(width) + x;
        zbuffer[loc] = z;
      }

      //bool l = false;

      void setColor(std::vector<Triangle>::iterator i, double x, double y, Point v1,
                        Point v2, Point v3)//int x, int y)
      {
        //This is probably where I want to do the interpolation, but the scanline maybe
        // cerr<<" ";
        //cerr<<i->getRed(0,0)<<i->getBlue(0,0)<<i->getGreen(0,0);

        if (x == this->width)
          {
            return;
          }
        else
        {


              //This might end up getting the wrong color you know, I think it will, I'm not sure why it hasn't been

          // cerr<<"SETTING COLOR"<<endl;
        //I think interpolating is basically slope
            //debugOt(v1,v2,v3,i);
          //Interpolate red
        double v5x = findx(v1,v2,y);
        // cerr<<"V5x "<<v5x<<endl;
        double v6x = findx(v2,v3,y);
        // cerr<<"V6x "<<v6x<<endl;


        //HAve to interpolate to get the Z value)
        double v5z = findc(v1.Z,v2.Z,v1.Y,v2.Y,y);
        double v6z = findc(v2.Z,v3.Z,v2.Y,v3.Y,y);
        double v4z = findc(v5z,v6z,v5x,v6x,x);

         if ((v4z <= getZ(x,y) ))//Not sure what to do when they ==
         {
         // cerr<<"DOES THIS HAPPEN"<<endl;
            return;//because if v4z <getZ(x,y) I don't actually want to draw the pixels
         }

//Guess I just need a semi-major refactor(again) to get the colors right

        // if (v4z < -1) cerr<<"V4Z is "<<v4z<<endl;
        // if (v1.Z < -1) cerr<<"V1Z is "<<v1.Z<<endl;
        // if (v2.Z < -1) cerr<<"V2Z is "<<v2.Z<<endl;
        // if (v3.Z < -1) cerr<<"V3Z is "<<v3.Z<<endl;
        setZ(x,y,v4z);

        double v5r = findc(v1.color[0],v2.color[0],v1.Y,v2.Y,y);
        double v6r = findc(v2.color[0],v3.color[0],v2.Y,v3.Y,y);
        double v4r = findc(v5r,v6r,v5x,v6x,x);
//that must mean it is interpolating the color to be 0 or something
        double v5u = findc(v1.color[2],v2.color[2],v1.Y,v2.Y,y);
        double v6u = findc(v2.color[2],v3.color[2],v2.Y,v3.Y,y);
        double v4u = findc(v5u,v6u,v5x,v6x,x);

        double v5g = findc(v1.color[1],v2.color[1],v1.Y,v2.Y,y);
        double v6g = findc(v2.color[1],v3.color[1],v2.Y,v3.Y,y);
        double v4g = findc(v5g,v6g,v5x,v6x,x);

           this->setRed(v4r, x, y);
           this->setBlue(v4u, x, y);
           this->setGreen(v4g, x, y);

//            if (x==332 && y == 321){
//         cerr<<"*****************"<<endl;
//         cerr<<"V1("<<v1.X<<","<<v1.Y<<","<<v1.Z<<")"<<endl;
//         cerr<<"V2("<<v2.X<<","<<v2.Y<<","<<v2.Z<<")"<<endl;
//         cerr<<"V3("<<v3.X<<","<<v3.Y<<","<<v3.Z<<")"<<endl;
//         cerr<<"v4("<<x<<","<<y<<")"<<endl;
//         cerr<<"v5x "<<v5x << " v6x "<<v6x<<endl;
//         cerr<<"v5r "<<v5r<<" v6r "<<v6r<<" v4r "<<endl;
//         cerr<<"v1r "<<i->getRed(0)<<" v2r "<<i->getRed(1)<<" v3r "<<i->getRed(2)<<endl;
// }
         // this->setRed(111,x,y);
         // this->setBlue(111,x,y);
         // this->setGreen(111,x,y);

        // cerr<<"V4R"<<v4r;


        }




        }
};

//Guy said somethings wrong with splitting of triangles and setting color to points, which sounds correct
std::vector<Triangle>
GetTriangles(void)
{
    vtkPolyDataReader *rdr = vtkPolyDataReader::New();
    rdr->SetFileName("proj1d_geometry.vtk");
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
    vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
    float *color_ptr = var->GetPointer(0);
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
        tris[idx].X[0] = pts->GetPoint(ptIds[0])[0];
        tris[idx].X[1] = pts->GetPoint(ptIds[1])[0];
        tris[idx].X[2] = pts->GetPoint(ptIds[2])[0];
        tris[idx].Y[0] = pts->GetPoint(ptIds[0])[1];
        tris[idx].Y[1] = pts->GetPoint(ptIds[1])[1];
        tris[idx].Y[2] = pts->GetPoint(ptIds[2])[1];
        tris[idx].Z[0] = pts->GetPoint(ptIds[0])[2];
        tris[idx].Z[1] = pts->GetPoint(ptIds[1])[2];
        tris[idx].Z[2] = pts->GetPoint(ptIds[2])[2];
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



void writeFlatBottom(Point v1, Point v2, Point v3, Screen screen,
              std::vector<Triangle>::iterator i)
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
            screen.setColor(i,ceil(leftend),y,v1,v2,v3);

        for (int x = ceil441(leftend); x <= floor441(rightend); x++)
          if (x >= 0 && x< screen.width)
            screen.setColor(i,x,y,v1,v2,v3);
      }
    }
}

void writeFlatTop(Point v1, Point v2, Point v3, Screen screen,
                    std::vector<Triangle>::iterator i)
{
    double leftslope, rightslope, li, ri,leftend,rightend;
    //debugOt(v1,v2,v3,i);
    /***Find the Left and Right Slope and Intercept***/

    //get rid of this
    //i->setBlue(255);
    //i->setGreen(0);
    //i->setRed(0);

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
          screen.setColor(i,ceil441(leftend),y,v1,v2,v3);

        for (int x = ceil441(leftend); x <= floor441(rightend); x++)
          if (x >= 0 && x< screen.width)
            screen.setColor(i,x,y,v1,v2,v3);
      }
    }
}


      void depositTriangle(std::vector<Triangle>::iterator i, Screen screen)
      {
        Point v1,v2,v3,v4;
        double slope, intercept;

        /*Step 1*/

    /*****I might need to do a continue if three vertices are in a line but maybe not*******/


        //3 Cases, when 1 and 2 are flat bottom and when 1 and 3 are flat bottom and when 2 and 3 are flat bottom
        if (i->getY1() == i->getY2())
        {
          //I think these can fuck up given a triangle like (0,0,1)(0,0,2)(5,2,0)
          //1 and 2 are the same Y
         // cerr<<"hi";
          if (i->getX1() < i->getX2())
          {
            i->setV(&v1,&v3,&v2);
          }
          else
          {
            i->setV(&v3,&v1,&v2);
          }
          if (v2.Y < v1.Y)
            writeFlatTop(v1,v2,v3,screen,i);
          else
            writeFlatBottom(v1,v2,v3,screen,i);
        }
        else if (i->getY2() == i->getY3())
        {
          //2 and 3 are the same Y
          if (i->getX2() < i->getX3())
          {
            i->setV(&v2,&v1,&v3);
          }
          else
          {
            //cerr<<"THIS SHOULD RUN";
            i->setV(&v2,&v3,&v1);
            //debugOt(v1,v2,v3,i);
          }
          if (v2.Y < v1.Y)
            writeFlatTop(v1,v2,v3,screen,i);
          else
            writeFlatBottom(v1,v2,v3,screen,i);
        }
        else if (i->getY1() == i->getY3())
        {
          //1 and 3 are the same Y,
          if (i->getX1() < i->getX3())
          {
            i->setV(&v1,&v2,&v3);
          }
          else
          {
            i->setV(&v3,&v2,&v1);
          }
          if (v2.Y < v1.Y)
            writeFlatTop(v1,v2,v3,screen,i);
          else
          {
            writeFlatBottom(v1,v2,v3,screen,i);
          }
        }
        else
        {
          //must split into a flat bottom and a flat top triangle

          /*Step 2*/


          i->sortVertices(&v1,&v2,&v3);
//cerr<<"V1"<<v1.X<<endl;
         // debugOt(v1,v2,v3,i);

          /*Step 3*/

          if (v1.X == v3.X) slope = 0; //right triangle
          else slope = (v1.Y - v3.Y)/(v1.X - v3.X);
          intercept = v3.Y - v3.X*slope;

          /*Step 4*/

          v4.Y = v2.Y;

          if (slope == 0) v4.X = v3.X;
          else v4.X = (v4.Y - intercept) / slope;

          /*Step 4.5*/
          //Gotta do some math to determine the Z for this point
          //I think I am just interpolating here


          //Wait isn't this going to be a different line based on where the split occurs???Probably some more mathy work
          //This isn't the correct way to interpolate but it is something.


          //How do I get a color for v4??????????????

          //If this doesn't work I need more if statements to determine exactly where point 4 is.

          //I guess the side the triangle is facing doesn't matter when interpolating the color
          v4.Z = findc(v1.Z,v3.Z,v1.Y,v3.Y,v4.Y);

          //since it's not drawing some, i think that means colors are undefined
          //sometimes
          v4.color[0] = findc(v1.color[0],v3.color[0],v1.Y,v3.Y,v4.Y);
          v4.color[1] = findc(v1.color[1],v3.color[1],v1.Y,v3.Y,v4.Y);
          v4.color[2] = findc(v1.color[2],v3.color[2],v1.Y,v3.Y,v4.Y);
          /*Step 5*/

          if (v2.X < v4.X){
            writeFlatTop(v2,v3,v4,screen,i);
          }
          else{
            writeFlatTop(v4,v3,v2,screen,i);
          }

          /*Step 6*/

          if (v2.X<v4.X){
            writeFlatBottom(v2,v1,v4,screen,i);
          }
          else {
            writeFlatBottom(v4,v1,v2,screen,i);
          }

    }
      }

int main()
{


  //Can I do this by having the screen have another buffer that is 1000/1000 and just stores the data, rather than fitting it all in
  //data structure?
   vtkImageData *image = NewImage(1000, 1000);
   unsigned char *buffer =
     (unsigned char *) image->GetScalarPointer(0,0,0);

   int npixels = 1000*1000;
    double *zbuffer = new double[npixels];
   for (int i = 0; i < npixels*3 ; i++)
       buffer[i] = 0; //You can use loop unrolling so this is not 2 for loops
   for (int i = 0; i < npixels; i++)
      zbuffer[i] = -1; //valid depth values go from -1 (back) to 0 (front)

   std::vector<Triangle> triangles = GetTriangles();

   Screen screen;
   screen.buffer = buffer;
   screen.zbuffer = zbuffer;
   screen.width = 1000;
   screen.height = 1000;


   double p1x, p1y, p2x, p2y, p3x, p3y;

   // YOUR CODE GOES HERE TO DEPOSIT TRIANGLES INTO PIXELS USING THE SCANLINE ALGORITHM
   int e = 0;
  for(std::vector<Triangle>::iterator i = triangles.begin(); i != triangles.end(); ++i)
  {
    //   if (e == 1) continue;

    // p1x = 50.0;
    // p1y = 90.0;
    // p2x = 40.0;
    // p2y = 40.0;
    // p3x = 60.0;
    // p3y = 40.0;
    //  i->setRed(1,0); i->setRed(0,1); i->setRed(0,2);
    //  i->setGreen(0,0); i->setGreen(1,1); i->setGreen(0,2);
    //  i->setBlue(0,0); i->setBlue(0,1); i->setBlue(1,2);
    //                                   //The depth is supposed to be between zero and 1, so something there is wonky??
    // i->setPT1(100,100,-.1);
    // i->setPT2(350,800,-.1);
    // i->setPT3(700,100,-.1);
    depositTriangle(i,screen);//I don't think I am assigning vertice 4 a color, that's actually something I know is happening
                                //and Is definitely the problem.

    //      i->setRed(0,0); i->setRed(0,1); i->setRed(0,2);
    //  i->setGreen(1,0); i->setGreen(1,1); i->setGreen(1,2);
    //  i->setBlue(0,0); i->setBlue(0,1); i->setBlue(0,2);

    // i->setPT1(200,200,0.5);
    // i->setPT2(500,200,0.5);
    // i->setPT3(400,300,-0.5);
    // depositTriangle(i,screen);

    //      i->setRed(1,0); i->setRed(1,1); i->setRed(1,2);
    //  i->setGreen(0,0); i->setGreen(0,1); i->setGreen(0,2);
    //  i->setBlue(0,0); i->setBlue(0,1); i->setBlue(0,2);

    //  i->setPT1(250,250,0.3);
    //  i->setPT2(350,350,0.3);
    //  i->setPT3(700,100,0.3);
    //  depositTriangle(i,screen);

    e++;
  }
  std::cerr<<e<<" ";


  std::cerr<<"Now IT's time To write";
   WriteImage(image, "allTriangles");
}
