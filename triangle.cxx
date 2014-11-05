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