//                                SPRINGDAMPER.CC                             


#include <springdamper.h>
#include <output.h>

s_springdamper::s_springdamper(void): C(0),K(0),Force(0),Length(0),
                                      Mode (_TRANSLATION),
                                      Ct(0),Kt(0),Torque(0),ANgle(0,RADIANS),
                                      Node1(0),Node2(0),
                                      _Node1(N),_Node2(N),_Mode(N),
                                      _C(N),_K(N),_Force(N),_Length(N),
                                      _Ct(N),_Kt(N),_Torque(N),_ANgle(N)
                                      {}

s_springdamper::~s_springdamper (void){}

inline const char* const s_springdamper::Gettype (void) const
{
   return "SPRINGDAMPER";
}

Boolean s_springdamper::Test(void)
{
   const int err_before=nerr;
   /* Controlla che siano presenti entrambi i marker */
   if (_Node1==N) out_error (34,"I");
   if (_Node2==N) out_error (34,"J");
   /* Controlla che non siano mischiati i tipi di definizione */
   if (Mode==_TRANSLATION) {
      if (_Ct==Y) out_error (35,"CT");
      if (_Kt==Y) out_error (35,"KT");
      if (_Torque==Y) out_error (35,"TORQUE");
      if (_ANgle==Y) out_error (35,"ANGLE");
   }
   if (Mode==_ROTATION) {
      if (_C==Y) out_error (36,"C");
      if (_K==Y) out_error (36,"K");
      if (_Force==Y) out_error (36,"FORCE");
      if (_Length==Y) out_error (36,"LENGTH");
   }
   if (err_before != nerr) return Y; else return N;
}

ostream& s_springdamper::Print (ostream& out) const
{
   out << endl 
       << "SPRINGDAMPER :" << label << "     TYPE:" << Mode << endl
       << "     Node 1 [" << _Node1 << "] = " << Node1 << endl
       << "     Node 2 [" << _Node2 << "] = " << Node2 << endl
       << "     C [" << _C << "] = " << C << endl
       << "     K [" << _K << "] = " << K << endl
       << "     FORCE [" << _Force << "] = " << Force << endl
       << "     LENGTH [" << _Length << "] = " << Length << endl
       << "     CT [" << _Ct << "] = " << Ct << endl
       << "     KT [" << _Kt << "] = " << Kt << endl
       << "     TORQUE [" << _Torque << "] = " << Torque << endl
       << "     ANGLE [" << _ANgle << "] = " << ANgle << endl;
   out << endl;
   return out;
}

void s_springdamper::Translate (ostream& out)
{
   return;
}
