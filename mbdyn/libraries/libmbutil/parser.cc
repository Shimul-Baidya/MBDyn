/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
 *
 * Pierangelo Masarati	<masarati@aero.polimi.it>
 * Paolo Mantegazza	<mantegazza@aero.polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 * Changing this copyright notice is forbidden.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/* parser */

#include <mbconfig.h>

#include <parser.h>

/* LowParser - begin */

static char
skip_remarks(InputStream& In)
{
skip_again:
   
   char cIn;
   while (isspace(cIn = In.get())) {
      NO_OP;
   }
   
   if (cIn == REMARK) {
      do {
	 NO_OP;
      } while ((cIn = In.get()) != '\n');
      goto skip_again;
   }
   
   if (cIn == '/') {
      if ((cIn = In.get()) == '*') {
	 for (;; cIn = In.get()) {
	    if (cIn == '*' && (cIn = In.get()) == '/') {	     
	       goto skip_again;
	    }	    
	 }	 
      } else {
	 In.putback(cIn);
	 return '/';
      }
   }
   
   return cIn;
}

void 
LowParser::PackWords(InputStream& In)
{
   char* pCur = sCurrWordBuf;
   char cIn = '\0';
   
   /* note: no remarks allowed inside words */
   while ((cIn = In.get()) != ':' && cIn != ';' && cIn != ',') {      
      if (!isspace(cIn)) {
	 *pCur++ = cIn;
	 if (pCur == sCurrWordBuf+iBufSize-1) {
	    *pCur = '\0';
	    return;
	 }
      }
   }
   
   *pCur = '\0';
   In.putback(cIn);
}


LowParser::Token 
LowParser::GetToken(InputStream& In)
{
   /* toglie gli spazi iniziali e tutti i commenti */
   char cIn = skip_remarks(In);
      
   if (isalpha(cIn) || cIn == '_') {
      PackWords(In.putback(cIn));
      return CurrToken = LowParser::WORD;
   }
   
   switch (cIn) {	     
    case ',':
      return CurrToken = LowParser::COMMA;
      
    case ':':
      return CurrToken = LowParser::COLON;       
      
    case ';':
      return CurrToken = LowParser::SEMICOLON;
      
    case '0':
    case '1':
    case '2':
    case '3':
    case '4':
    case '5':
    case '6':
    case '7':
    case '8':
    case '9':
    case '.':
    case '-':
    case '+':
      In.putback(cIn) >> dCurrNumber;
      return CurrToken = LowParser::NUMBER;
      
    default:
      In.putback(cIn);
      return CurrToken = LowParser::UNKNOWN;      
   }
}


doublereal 
LowParser::dGetReal(void)
{
   return dCurrNumber;
}


integer 
LowParser::iGetInt(void)
{
   return integer(dCurrNumber);
}


char* 
LowParser::sGetWord(void)
{
   return sCurrWordBuf;
}   

/* LowParser - end */


/* KeyTable - begin */

KeyTable::KeyTable(int iTableLen, const char* const sTable[])
{
   sKeyWords = (char* const*)sTable;
   iNumKeys = iTableLen;
}
 

int 
KeyTable::Find(const char* sToFind)
{
   int iCnt = -1;
   while (iCnt++,  iCnt < iNumKeys) {      
      if (strcasecmp(sKeyWords[iCnt], sToFind) == 0) {	 
	 return iCnt;
      }
   }   
   
   return -1;
}  

/* KeyTable - end */


/* HighParser - begin */

HighParser::HighParser(MathParser& MP, KeyTable& KT, InputStream& streamIn)
: ESCAPE_CHAR('\\'),
pIn(&streamIn), pf(NULL),
MathP(MP), 
KeyT(KT)
{
   DEBUGCOUTFNAME("HighParser::HighParser");
   CurrToken = HighParser::DESCRIPTION;  
}   


HighParser::~HighParser(void)
{
   DEBUGCOUTFNAME("HighParser::~HighParser");
   Close();
}
 

void 
HighParser::Close(void)
{
   NO_OP;
}


void 
HighParser::PutKeyTable(KeyTable& KT)
{
   KeyT = KT;
}

MathParser& 
HighParser::GetMathParser(void) 
{
   return MathP;
}

int 
HighParser::GetLineNumber(void) const
{
   return ((InputStream*)pIn)->GetLineNumber();
}


HighParser::ErrOut 
HighParser::GetLineData(void) const
{
   ErrOut LineData;
   LineData.iLineNumber = GetLineNumber();
   LineData.sFileName = NULL;
   return LineData;	
}


flag 
HighParser::fIsDescription(void)
{
   if (CurrToken != HighParser::DESCRIPTION) {
      cerr << "Parser error in HighParser::fIsDescription, invalid call to GetDescription at line "
	<< GetLineData() << endl;
      return flag(0);
   }
   return flag(1);
}


int 
HighParser::iGetDescription_(const char* const s)
{
   int i = KeyT.Find(s);
   
   if (FirstToken() == HighParser::UNKNOWN) {
      cerr << "Parser error in HighParser::iGetDescription_(), semicolon expected at line " 
	<< GetLineData() << endl;
      THROW(HighParser::ErrSemicolonExpected());      
   }
   
   return i;   
}


void 
HighParser::Set_(void) 
{      
   if (FirstToken() == UNKNOWN) {
      cerr << endl << "Parser error in HighParser::Set_(), colon expected at line " 
	<< GetLineData() << endl;
      THROW(HighParser::ErrColonExpected());
   }
   
   GetReal();   
}


int 
HighParser::GetDescription(void)
{
   const char sFuncName[] = "HighParser::GetDescription()";
   
   /* Checks if current token is a description */
   if (!fIsDescription()) {
      THROW(HighParser::ErrInvalidCallToGetDescription());
   }      
   
restart_parsing:
   
   if ((CurrLowToken = LowP.GetToken(*pIn)) != LowParser::WORD) {
      if (pIn->GetStream().eof()) {
	 THROW(ErrFile());
      } else {     	 
	 cerr << endl << "Parser error in "
	   << sFuncName << ", keyword expected at line " 
	   << GetLineData() << endl;
	 THROW(HighParser::ErrKeyWordExpected());
      }      
   }
   
   /* Description corrente */
   char* s = LowP.sGetWord();
   
   if (strcmp(s, "set") == 0) {
      Set_();
      goto restart_parsing;
   } /* else */
   return iGetDescription_(s);
}


HighParser::Token 
HighParser::FirstToken(void)
{
   if ((CurrLowToken = LowP.GetToken(*pIn)) == LowParser::COLON) {
      return (CurrToken = HighParser::ARG);
   } else {
      if (CurrLowToken != LowParser::SEMICOLON) {
	 return (CurrToken = HighParser::UNKNOWN);
      } else {
	 return (CurrToken = HighParser::DESCRIPTION);
      }
   }
}

void 
HighParser::ExpectDescription(void)
{
   CurrToken = HighParser::DESCRIPTION;
}


void 
HighParser::ExpectArg(void)
{
   CurrToken = HighParser::ARG;
}


flag 
HighParser::fIsArg(void)
{
   if (CurrToken == ARG) {
      return flag(1);
   }   
   return flag(0);
}

void 
HighParser::PutBackSemicolon(void)
{
   if (CurrLowToken == LowParser::SEMICOLON) {
      pIn->putback(';');
   }   
}


void 
HighParser::NextToken(const char* sFuncName)
{
   switch (CurrLowToken = LowP.GetToken(*pIn)) {
    case LowParser::COMMA:
      CurrToken = HighParser::ARG;
      break;
    case LowParser::SEMICOLON:
      CurrToken = HighParser::DESCRIPTION;
      break;
    default:
      cerr << endl << "Parser error in "
	<< sFuncName << ", missing separator at line " 
	<< GetLineData() << endl;
      THROW(HighParser::ErrMissingSeparator());
   }   
}


flag 
HighParser::IsKeyWord(const char* sKeyWord)
{      
   const char sFuncName[] = "HighParser::IsKeyWord()";
      
   if (CurrToken != HighParser::ARG) {      
      return 0;
   }   
   
   char* sBuf = sStringBuf;
   char* sBufWithSpaces = sStringBufWithSpaces;

   char cIn = skip_remarks(*pIn);
   
   if (!isalpha(cIn)) {
      pIn->putback(cIn);
      return 0;
   }
   
   *sBuf++ = cIn;
   *sBufWithSpaces++ = cIn;
   
   /* Forse e' meglio modificare in modo che digerisca anche gli spazi
    * tra due parole, magari con due buffer, uno in cui li mangia per fare
    * il confronto con la keyword, l'altro in cui li tiene per l'eventuale 
    * putback */
   while (isalnum(cIn = pIn->get()) || isspace(cIn)) {
      *sBufWithSpaces++ = cIn;
      if (isalnum(cIn)) {
	 *sBuf++ = cIn;
      }
      if (sBufWithSpaces >= sStringBufWithSpaces+iBufSize-1) {
	 break;
      }      
   }   
   pIn->putback(cIn);
   
   *sBuf = '\0';
   *sBufWithSpaces = '\0';
   
   if (!strcasecmp(sStringBuf, sKeyWord)) {
      NextToken(sFuncName);		
      return 1;
   }   
   
   while (sBufWithSpaces > sStringBufWithSpaces) {      
      pIn->putback(*--sBufWithSpaces);
   }   
   
   return 0;
}


int 
HighParser::IsKeyWord(void)
{
   const char sFuncName[] = "HighParser::IsKeyWord()";
   
   if (CurrToken != HighParser::ARG) {      
      return -1;
   }
      
   char* sBuf = sStringBuf;
   char* sBufWithSpaces = sStringBufWithSpaces;

   char cIn = skip_remarks(*pIn);
   
   if (!isalpha(cIn)) {
      pIn->putback(cIn);
      return -1;
   }   
   *sBuf++ = cIn;
   *sBufWithSpaces++ = cIn;
   
   while (isalnum(cIn = pIn->get()) || isspace(cIn)) {
      *sBufWithSpaces++ = cIn;
      if (isalnum(cIn)) {
	 *sBuf++ = cIn;
      }
      if (sBufWithSpaces >= sStringBufWithSpaces+iBufSize-1) {
	 break;
      }      
   }   
   pIn->putback(cIn);
   
   *sBuf = '\0';
   *sBufWithSpaces = '\0';
   
   int iKW;
   if ((iKW = KeyT.Find(sStringBuf)) >= 0) {
      NextToken(sFuncName);
      return iKW;
   }   
   
   while (sBufWithSpaces > sStringBufWithSpaces) {
      pIn->putback(*--sBufWithSpaces);
   }   
   
   return -1;   
}


integer 
HighParser::GetInt(int iDefval)
{
   const char sFuncName[] = "HighParser::GetInt()";
   
   if (CurrToken != HighParser::ARG) {
      cerr << endl << "Parser error in "
	<< sFuncName << ", integer arg expected at line " 
	<< GetLineData() << endl;
      THROW(HighParser::ErrIntegerExpected());
   }
   
   integer iReturnValue;
   
#ifdef USE_EXCEPTIONS
   try {
#endif
      
      iReturnValue = int(MathP.Get(*pIn, (double)iDefval));
      
#ifdef USE_EXCEPTIONS
   }
   catch (MathParser::ErrGeneric e) {
      cerr << sFuncName << ": error return from MathParser at line " << GetLineData() << endl;
      throw e;
   }
#endif     
   
   NextToken(sFuncName);
   return iReturnValue;
}


doublereal 
HighParser::GetReal(double dDefval)
{
   const char sFuncName[] = "HighParser::GetReal()";
   
   if (CurrToken != HighParser::ARG) {
      cerr << endl << "Parser error in "
	<< sFuncName << ", real arg expected at line " 
	<< GetLineData() << endl;
      THROW(HighParser::ErrRealExpected());
   }
   
   doublereal dReturnValue;
#ifdef USE_EXCEPTIONS
   try {
#endif
   
      dReturnValue = doublereal(MathP.Get(*pIn, dDefval));

#ifdef USE_EXCEPTIONS
   }
   catch (MathParser::ErrGeneric e) {
      cerr << sFuncName << ": error return from MathParser at line " << GetLineData() << endl;
      throw e;
   }
#endif     
      
   NextToken(sFuncName);
   return dReturnValue;
}


int 
HighParser::GetWord(void)
{
   const char sFuncName[] = "HighParser::GetWord()";
   
   if (CurrToken != HighParser::ARG) {
      cerr << endl << "Parser error in "
	<< sFuncName << ", keyword arg expected at line " 
	<< GetLineData() << endl;
      THROW(HighParser::ErrKeyWordExpected());
   }
   
   if ((CurrLowToken = LowP.GetToken(*pIn)) != LowParser::WORD) {
      cerr << endl << "Parser error in "
	<< sFuncName << ", keyword expected at line " 
	<< GetLineData() << endl;
      THROW(HighParser::ErrKeyWordExpected());
   }
   
   int i = KeyT.Find(LowP.sGetWord());
   
   NextToken(sFuncName);
   return i;
}


const char* 
HighParser::GetString(void)
{
   const char sFuncName[] = "HighParser::GetString()";
   
   cerr << "line " << GetLineData()
     << ": warning, use of deprecated method \"GetString\"" << endl;
   
   if (CurrToken != HighParser::ARG) {
      cerr << endl << "Parser error in "
	<< sFuncName << ", string arg expected at line " 
	<< GetLineData() << endl;
      THROW(HighParser::ErrStringExpected());
   }
   
   char* s = sStringBuf;
   char* sTmp = s;
   
   char cIn = '\0';
   
   while (isspace(cIn = pIn->get())) {
      NO_OP;
   }   
   
   pIn->putback(cIn);
   while ((cIn = pIn->get()) != ',' && cIn != ';') {
      /* Attenzione! cosi' la legge tutta, 
       * ma ne tiene solo iBufSize-1 caratteri */
      if (sTmp < s+iBufSize-1) {
	 *sTmp++ = cIn;
      }      
   }
   
   pIn->putback(cIn);
   *sTmp = '\0';	
   
   NextToken(sFuncName);
   return s;
}


const char* 
HighParser::GetStringWithDelims(enum Delims Del)
{
   const char sFuncName[] = "HighParser::GetStringWithDelims()";
   
   if (CurrToken != HighParser::ARG) {
      cerr << endl << "Parser error in "
	<< sFuncName << ", string arg expected at line " 
	<< GetLineData() << endl;
      THROW(HighParser::ErrStringExpected());
   }
   
   char* s = sStringBuf;
   char* sTmp = s;
   
   char cLdelim = '\0';
   char cRdelim = '\0';
   
   switch (Del) {
    case PLAINBRACKETS:
      cLdelim = '(';
      cRdelim = ')';
      break;
    case SQUAREBRACKETS:
      cLdelim = '[';
      cRdelim = ']';
      break;
    case CURLYBRACKETS:
      cLdelim = '{';
      cRdelim = '}';
      break;
    case SINGLEQUOTE:
      cLdelim = '`';
      cRdelim = '\'';
      break;
    default:
    case UNKNOWNDELIM:
    case DEFAULTDELIM:
    case DOUBLEQUOTE:
      cLdelim = '"';
      cRdelim = '"';
      break;	
   }

   char cIn = skip_remarks(*pIn);
      
   /* Se trova il delimitatore sinistro, legge la stringa */
   if (cIn == cLdelim) {
      while ((cIn = pIn->get()) != cRdelim) {
	 /* Attenzione! cosi' la legge tutta, 
	  * ma ne tiene solo iBufSize-1 caratteri */
	 if (sTmp < s+iBufSize-1) {
	    if (cIn == ESCAPE_CHAR) {
	       *sTmp++ = cIn;
	       cIn = pIn->get();
	    }
	    *sTmp++ = cIn;	   
	 }	 
      }   
   
   /* Se trova una virgola o un punto e virgola, le rimette nello stream e passa oltre, 
    * restituendo un puntatore nullo. Il chiamante deve occuparsi della
    * gestione del valore di default */
   } else if (cIn == ',' || cIn == ';') {	
      pIn->putback(cIn);
      goto nullstring;  
   
   /* Altrimenti c'e' qualcosa senza delimitatore. Adesso da' errore,
    * forse e' piu' corretto fargli ritornare lo stream intatto */
   } else {	
      cerr << endl << "Parser error in "
	<< sFuncName << endl 
	<< "first non-blank char at line " 
	<< GetLineData() << " isn't a valid left-delimiter;" << endl
	<< "aborting ..." << endl;
      THROW(HighParser::ErrIllegalDelimiter());
   }   
   
   /* Mette zero al termine della stringa */
   *sTmp = '\0';	
	
   nullstring:
   
   NextToken(sFuncName);
   return s;
}


/* Legge un Vec3 */
Vec3 
HighParser::GetVec3(void)
{
   Vec3 v(0.);
   return GetVec3(v);
}


/* Legge un Vec3 */
Vec3 
HighParser::GetVec3(const Vec3& vDef)
{
   if (IsKeyWord("null")) {	
      return vDef;
   }      
   
   doublereal x1 = GetReal(vDef.dGet(1));
   doublereal x2 = GetReal(vDef.dGet(2));
   doublereal x3 = GetReal(vDef.dGet(3));
   return Vec3(x1, x2, x3);
}


/* Legge una matrice R sotto forma di due vettori (oppure eye) */
Mat3x3 
HighParser::GetMatR2vec(void)
{
   if (IsKeyWord("eye")) {	
      return Eye3;
   }
   
   if (IsKeyWord("matr")) {
      doublereal r11 = GetReal();
      doublereal r12 = GetReal();
      doublereal r13 = GetReal();
      doublereal r21 = GetReal();
      doublereal r22 = GetReal();
      doublereal r23 = GetReal();
      doublereal r31 = GetReal();
      doublereal r32 = GetReal();
      doublereal r33 = GetReal();
      
      return Mat3x3(r11, r21, r31, r12, r22, r32, r13, r23, r33);
   }
   
   if (IsKeyWord("euler")) {
      doublereal e1 = GetReal();
      doublereal e2 = GetReal();
      doublereal e3 = GetReal();
      
      return RFromEulerAngles(Vec3(e1, e2, e3));
   }
   
   int i1 = GetInt();
   doublereal x1 = GetReal();
   doublereal x2 = GetReal();
   doublereal x3 = GetReal();
   Vec3 v1(x1, x2, x3);
   
   int i2 = GetInt();
   x1 = GetReal();
   x2 = GetReal();
   x3 = GetReal();
   Vec3 v2(x1, x2, x3);    
   
   return MatR2vec(i1, v1, i2, v2);
}


/* Legge una matrice 3x3 simmetrica come diagonale o triangolare superiore */
Mat3x3 
HighParser::GetMat3x3Sym(void)
{
   if (IsKeyWord("null")) {	
      return Mat3x3(0.);
   } else if (IsKeyWord("eye")) {	
      return Eye3;
   } else if (IsKeyWord("diag")) {	
      doublereal m11 = GetReal();
      doublereal m22 = GetReal();
      doublereal m33 = GetReal();
      return Mat3x3(m11, 0., 0., 0., m22, 0., 0., 0., m33);
   } else {
      doublereal m11 = GetReal();
      doublereal m12 = GetReal();
      doublereal m13 = GetReal();
      doublereal m22 = GetReal();
      doublereal m23 = GetReal();
      doublereal m33 = GetReal();
      return Mat3x3(m11, m12, m13, m12, m22, m23, m13, m23, m33);
   }   
}

/* Legge una matrice 3x3 generica (diagonale o nulla) */
Mat3x3 
HighParser::GetMat3x3(void)
{
   Mat3x3 m(0.);
   return GetMat3x3(m);
}


/* Legge una matrice 3x3 generica (diagonale o nulla) */
Mat3x3 
HighParser::GetMat3x3(const Mat3x3& mDef)
{
   if (IsKeyWord("null")) {	
      return mDef;
   } else if(IsKeyWord("eye")) {	
      return Eye3;
   } else if(IsKeyWord("diag")) {	
      doublereal m11 = GetReal(mDef.dGet(1, 1));
      doublereal m22 = GetReal(mDef.dGet(2, 2));
      doublereal m33 = GetReal(mDef.dGet(3, 3));
      return Mat3x3(m11, 0., 0., 0., m22, 0., 0., 0., m33);
   } else if(IsKeyWord("sym")) {	
      doublereal m11 = GetReal(mDef.dGet(1, 1));
      doublereal m12 = GetReal(mDef.dGet(1, 2));
      doublereal m13 = GetReal(mDef.dGet(1, 3));
      doublereal m22 = GetReal(mDef.dGet(2, 2));
      doublereal m23 = GetReal(mDef.dGet(2, 3));
      doublereal m33 = GetReal(mDef.dGet(3, 3));
      return Mat3x3(m11, m12, m13, m12, m22, m23, m13, m23, m33);
   } else {
      doublereal m11 = GetReal(mDef.dGet(1, 1));
      doublereal m12 = GetReal(mDef.dGet(1, 2));
      doublereal m13 = GetReal(mDef.dGet(1, 3));
      doublereal m21 = GetReal(mDef.dGet(2, 1));
      doublereal m22 = GetReal(mDef.dGet(2, 2));
      doublereal m23 = GetReal(mDef.dGet(2, 3));
      doublereal m31 = GetReal(mDef.dGet(3, 1));
      doublereal m32 = GetReal(mDef.dGet(3, 2));
      doublereal m33 = GetReal(mDef.dGet(3, 3));
      return Mat3x3(m11, m21, m31, m12, m22, m32, m13, m23, m33);
   }   
}


/* Legge un Vec6 */
Vec6 
HighParser::GetVec6(void)
{
   Vec6 v(0.);
   return GetVec6(v);
}


/* Legge un Vec6 */
Vec6 
HighParser::GetVec6(const Vec6& vDef)
{
   if (IsKeyWord("null")) {	
      return vDef;
   }      
   
   doublereal x1 = GetReal(vDef.dGet(1));
   doublereal x2 = GetReal(vDef.dGet(2));
   doublereal x3 = GetReal(vDef.dGet(3));
   doublereal x4 = GetReal(vDef.dGet(4));
   doublereal x5 = GetReal(vDef.dGet(5));
   doublereal x6 = GetReal(vDef.dGet(6));
   return Vec6(x1, x2, x3, x4, x5, x6);
}


/* Legge una matrice 6x6 generica (diagonale o nulla) */
Mat6x6 
HighParser::GetMat6x6(void)
{
   Mat6x6 m(0.);
   return GetMat6x6(m);
}


/* Legge una matrice 6x6 generica (diagonale o nulla) */
Mat6x6 
HighParser::GetMat6x6(const Mat6x6& mDef)
{
   if (IsKeyWord("null")) {	
      return Zero6x6;
   } else if(IsKeyWord("eye")) {	
      return Eye6;
   } else if(IsKeyWord("diag")) {	
      doublereal m11 = GetReal(mDef.dGet(1, 1));
      doublereal m22 = GetReal(mDef.dGet(2, 2));
      doublereal m33 = GetReal(mDef.dGet(3, 3));
      doublereal m44 = GetReal(mDef.dGet(4, 4));
      doublereal m55 = GetReal(mDef.dGet(5, 5));
      doublereal m66 = GetReal(mDef.dGet(6, 6));
      return Mat6x6(m11, 0., 0., 0., 0., 0., 
		    0., m22, 0., 0., 0., 0.,
		    0., 0., m33, 0., 0., 0.,
		    0., 0., 0., m44, 0., 0.,
		    0., 0., 0., 0., m55, 0.,
		    0., 0., 0., 0., 0., m66);
   } else if(IsKeyWord("sym")) {	
      doublereal m11 = GetReal(mDef.dGet(1, 1));
      doublereal m12 = GetReal(mDef.dGet(1, 2));
      doublereal m13 = GetReal(mDef.dGet(1, 3));
      doublereal m14 = GetReal(mDef.dGet(1, 4));
      doublereal m15 = GetReal(mDef.dGet(1, 5));
      doublereal m16 = GetReal(mDef.dGet(1, 6));
      doublereal m22 = GetReal(mDef.dGet(2, 2));
      doublereal m23 = GetReal(mDef.dGet(2, 3));
      doublereal m24 = GetReal(mDef.dGet(2, 4));
      doublereal m25 = GetReal(mDef.dGet(2, 5));
      doublereal m26 = GetReal(mDef.dGet(2, 6));
      doublereal m33 = GetReal(mDef.dGet(3, 3));
      doublereal m34 = GetReal(mDef.dGet(3, 4));
      doublereal m35 = GetReal(mDef.dGet(3, 5));
      doublereal m36 = GetReal(mDef.dGet(3, 6));
      doublereal m44 = GetReal(mDef.dGet(4, 4));
      doublereal m45 = GetReal(mDef.dGet(4, 5));
      doublereal m46 = GetReal(mDef.dGet(4, 6));
      doublereal m55 = GetReal(mDef.dGet(5, 5));
      doublereal m56 = GetReal(mDef.dGet(5, 6));
      doublereal m66 = GetReal(mDef.dGet(6, 6));      
      return Mat6x6(m11, m12, m13, m14, m15, m16, 
		    m12, m22, m23, m24, m25, m26,
		    m13, m23, m33, m34, m35, m36,
		    m14, m24, m34, m44, m45, m46,
		    m15, m25, m35, m45, m55, m56,
		    m16, m26, m36, m46, m56, m66);
   } else if(IsKeyWord("anba")) {
      /* Formato ANBA, in cui vale la trasformazione:
       * ex = e2
       * ey = e3
       * ez = e1
       */
      doublereal m22 = GetReal(mDef.dGet(2, 2));
      doublereal m23 = GetReal(mDef.dGet(2, 3));
      doublereal m21 = GetReal(mDef.dGet(2, 1));
      doublereal m25 = GetReal(mDef.dGet(2, 5));
      doublereal m26 = GetReal(mDef.dGet(2, 6));
      doublereal m24 = GetReal(mDef.dGet(2, 4));
      doublereal m32 = GetReal(mDef.dGet(3, 2));
      doublereal m33 = GetReal(mDef.dGet(3, 3));
      doublereal m31 = GetReal(mDef.dGet(3, 1));
      doublereal m35 = GetReal(mDef.dGet(3, 5));
      doublereal m36 = GetReal(mDef.dGet(3, 6));
      doublereal m34 = GetReal(mDef.dGet(3, 4));
      doublereal m12 = GetReal(mDef.dGet(1, 2));
      doublereal m13 = GetReal(mDef.dGet(1, 3));
      doublereal m11 = GetReal(mDef.dGet(1, 1));
      doublereal m15 = GetReal(mDef.dGet(1, 5));
      doublereal m16 = GetReal(mDef.dGet(1, 6));
      doublereal m14 = GetReal(mDef.dGet(1, 4));
      doublereal m52 = GetReal(mDef.dGet(5, 2));
      doublereal m53 = GetReal(mDef.dGet(5, 3));
      doublereal m51 = GetReal(mDef.dGet(5, 1));
      doublereal m55 = GetReal(mDef.dGet(5, 5));
      doublereal m56 = GetReal(mDef.dGet(5, 6));
      doublereal m54 = GetReal(mDef.dGet(5, 4));
      doublereal m62 = GetReal(mDef.dGet(6, 2));
      doublereal m63 = GetReal(mDef.dGet(6, 3));
      doublereal m61 = GetReal(mDef.dGet(6, 1));
      doublereal m65 = GetReal(mDef.dGet(6, 5));
      doublereal m66 = GetReal(mDef.dGet(6, 6));
      doublereal m64 = GetReal(mDef.dGet(6, 4));
      doublereal m42 = GetReal(mDef.dGet(4, 2));
      doublereal m43 = GetReal(mDef.dGet(4, 3));
      doublereal m41 = GetReal(mDef.dGet(4, 1));
      doublereal m45 = GetReal(mDef.dGet(4, 5));
      doublereal m46 = GetReal(mDef.dGet(4, 6));
      doublereal m44 = GetReal(mDef.dGet(4, 4));

      return Mat6x6(m11, m21, m31, m41, m51, m61,
		    m12, m22, m32, m42, m52, m62,
		    m13, m23, m33, m43, m53, m63,
		    m14, m24, m34, m44, m54, m64,
		    m15, m25, m35, m45, m55, m65,
		    m16, m26, m36, m46, m56, m66);
   } else {
      doublereal m11 = GetReal(mDef.dGet(1, 1));
      doublereal m12 = GetReal(mDef.dGet(1, 2));
      doublereal m13 = GetReal(mDef.dGet(1, 3));
      doublereal m14 = GetReal(mDef.dGet(1, 4));
      doublereal m15 = GetReal(mDef.dGet(1, 5));
      doublereal m16 = GetReal(mDef.dGet(1, 6));
      doublereal m21 = GetReal(mDef.dGet(2, 1));
      doublereal m22 = GetReal(mDef.dGet(2, 2));
      doublereal m23 = GetReal(mDef.dGet(2, 3));
      doublereal m24 = GetReal(mDef.dGet(2, 4));
      doublereal m25 = GetReal(mDef.dGet(2, 5));
      doublereal m26 = GetReal(mDef.dGet(2, 6));
      doublereal m31 = GetReal(mDef.dGet(3, 1));
      doublereal m32 = GetReal(mDef.dGet(3, 2));
      doublereal m33 = GetReal(mDef.dGet(3, 3));
      doublereal m34 = GetReal(mDef.dGet(3, 4));
      doublereal m35 = GetReal(mDef.dGet(3, 5));
      doublereal m36 = GetReal(mDef.dGet(3, 6));
      doublereal m41 = GetReal(mDef.dGet(4, 1));
      doublereal m42 = GetReal(mDef.dGet(4, 2));
      doublereal m43 = GetReal(mDef.dGet(4, 3));
      doublereal m44 = GetReal(mDef.dGet(4, 4));
      doublereal m45 = GetReal(mDef.dGet(4, 5));
      doublereal m46 = GetReal(mDef.dGet(4, 6));
      doublereal m51 = GetReal(mDef.dGet(5, 1));
      doublereal m52 = GetReal(mDef.dGet(5, 2));
      doublereal m53 = GetReal(mDef.dGet(5, 3));
      doublereal m54 = GetReal(mDef.dGet(5, 4));
      doublereal m55 = GetReal(mDef.dGet(5, 5));
      doublereal m56 = GetReal(mDef.dGet(5, 6));
      doublereal m61 = GetReal(mDef.dGet(6, 1));
      doublereal m62 = GetReal(mDef.dGet(6, 2));
      doublereal m63 = GetReal(mDef.dGet(6, 3));
      doublereal m64 = GetReal(mDef.dGet(6, 4));
      doublereal m65 = GetReal(mDef.dGet(6, 5));
      doublereal m66 = GetReal(mDef.dGet(6, 6));
      return Mat6x6(m11, m21, m31, m41, m51, m61,
		    m12, m22, m32, m42, m52, m62,
		    m13, m23, m33, m43, m53, m63,
		    m14, m24, m34, m44, m54, m64,
		    m15, m25, m35, m45, m55, m65,
		    m16, m26, m36, m46, m56, m66);
   }   
}

/* provvisoria */
void 
HighParser::GetMat6xN(Mat3xN& m1, Mat3xN& m2, integer iNumCols)
{
   ASSERT(iNumCols > 0);
   ASSERT(m1.iGetNumCols() == iNumCols);
   ASSERT(m2.iGetNumCols() == iNumCols);
   
   if(IsKeyWord("null")) {
      m1.Reset();
      m2.Reset();
   } else if (IsKeyWord("anba")) {
      int vi[] = { 2, 3, 1 };
      for (int i = 0; i < 3; i++) {
	 for (integer j = 1; j <= iNumCols; j++) {
	    m1.Put(vi[i], j, GetReal());
	 }
      }      
      for (int i = 0; i < 3; i++) {
	 for (integer j = 1; j <= iNumCols; j++) {
	    m2.Put(vi[i], j, GetReal());
	 }
      }      
   } else {     
      for (int i = 1; i <= 3; i++) {
	 for (integer j = 1; j <= iNumCols; j++) {
	    m1.Put(i, j, GetReal());
	 }
      }      
      for (int i = 1; i <= 3; i++) {
	 for (integer j = 1; j <= iNumCols; j++) {
	    m2.Put(i, j, GetReal());
	 }
      }
   }
}


/* HighParser - end */

ostream& 
operator << (ostream& out, const HighParser::ErrOut& err)
{
   out << err.iLineNumber;
   if (err.sFileName != NULL) {      
      out << ", file <" << err.sFileName << '>';
   }
   return out;   
}
