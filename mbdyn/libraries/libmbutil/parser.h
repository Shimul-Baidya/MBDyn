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
 * the Free Software Foundation (version 2 of the License).
 * 
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

/* Parser per l'ingresso dati - parte generale */

/* Si compone di tre diverse strutture di scansione, 
 * piu' le strutture di memorizzazione di stringhe e variabili.
 * 
 * La prima struttura e' data dal LowParser, che riconosce gli elementi 
 * della sintassi. E' data da:
 * <statement_list>::=
 *   <statement> ; <statement_list>
 *   epsilon
 * <statement>::=
 *   <description>
 *   <description> : <arg_list>
 * <arg_list>::=
 *   <arg>
 *   <arg> , <arg_list>
 * <arg>::=
 *   <word>
 *   <number>
 * ecc. ecc.
 * 
 * La seconda struttura e' data dall'HighParser, che riconosce la sintassi 
 * vera e propria. In alternativa al LowParser, qualora sia atteso un valore
 * numerico esprimibile mediante un'espressione regolare 
 * (espressione matematica), e' possibile invocare la terza struttura,
 * il MathParser. Questo analizza espressioni anche complesse e multiple,
 * se correttamente racchiuse tra parentesi.
 * 
 * L'HighParser deve necessariamente riconoscere una parola chiave nel campo
 * <description>, mentre puo' trovare parole qualsiasi nel campo <arg>
 * qualora sia attesa una stringa.
 * 
 * Le parole chiave vengono fornite all'HighParser attraverso la KeyTable, 
 * ovvero una lista di parole impaccate (senza spazi). L'uso consigliato e':
 * 
 *   const char sKeyWords[] = { "keyword0",
 *                              "keyword1",
 *                              "...",
 *                              "keywordN"};
 * 
 *   enum KeyWords { KEYWORD0 = 0,
 *                   KEYWORD1,
 *                   ...,
 *                   KEYWORDN,
 *                   LASTKEYWORD};
 * 
 *   KeyTable K((int)LASTKEYWORD, sKeyWords);
 * 
 * Il MathParser usa una tabella di simboli, ovvero nomi (dotati di tipo) 
 * a cui e' associato un valore. La tabella e' esterna e quindi puo' essere 
 * conservata ed utilizzata in seguito conservando in memoria i nomi
 * definiti in precedenza.
 * 
 * A questo punto si puo' generare la tabella dei simboli:
 * 
 *   int iSymbolTableInitialSize = 10;
 *   Table T(iSymbolTableInitialSize);
 * 
 * Quindi si crea il MathParser:
 * 
 *   MathParser Math(T);
 * 
 * Infine si genera l'HighParser:
 * 
 *   HighParser HP(Math, K, StreamIn);
 * 
 * dove StreamIn e' l'istream da cui avviene la lettura.
 */
            

#ifndef PARSER_H
#define PARSER_H


#include <iostream.h>
#include <fstream.h>
extern "C" {
#include <strings.h>
#include <ctype.h>
#include <myf2c.h>

#include <stdlib.h>
#include <unistd.h>
}
   
#include <input.h>
#include <mathp.h>
#include <matvec3.h>
#include <matvec3n.h>
#include <matvec6.h>


/* Classi dichiarate */
class LowParser;
class KeyTable;
class HighParser;


const unsigned int iBufSize = 1024;


/* LowParser - begin */

class LowParser {

   friend class HighParser;

 public:
   enum Token {
      UNKNOWN,
	WORD,
	COMMA = ',',
	COLON = ':',
	SEMICOLON = ';',
	NUMBER
   };
   
 private:
   enum Token CurrToken;
   char sCurrWordBuf[iBufSize];
   doublereal dCurrNumber;

   void PackWords(InputStream& In);
   
 public:
   Token GetToken(InputStream& In);      
   doublereal dGetReal(void);
   integer iGetInt(void);   
   char* sGetWord(void);      
};

/* LowParser - end */


/* KeyTable - begin */

class KeyTable {
 private:
   char* const* sKeyWords;
   int iNumKeys;
   
 public:      
   KeyTable(int iTableLen, const char* const sTable[]);   
   int Find(const char* sToFind);  
   
};

/* KeyTable - end */


/* HighParser - begin */

class HighParser {

 public:   
   class ErrInvalidCallToGetDescription {};
   class ErrKeyWordExpected {};
   class ErrSemicolonExpected {};
   class ErrColonExpected {};
   class ErrMissingSeparator {};
   class ErrIntegerExpected {};
   class ErrRealExpected {};
   class ErrStringExpected {};
   class ErrIllegalDelimiter {};
   
 public:      
   enum Token {
      UNKNOWN = -1,
	DESCRIPTION,
	FIRSTARG,
	ARG,
	LASTARG,
	NOARGS,
	WORD,	
	NUMBER,
	STRING,
	LASTITEM
   };
   
   enum Delims {
      UNKNOWNDELIM = -1,
	PLAINBRACKETS,
	SQUAREBRACKETS,
	CURLYBRACKETS,
	SINGLEQUOTE,
	DOUBLEQUOTE,
	DEFAULTDELIM,
	LASTDELIM
   };
   
   const char ESCAPE_CHAR;
   
 public:
   struct ErrOut {
      const char* sFileName;
      unsigned int iLineNumber;
   };
   
 protected:   
   /* Parser di basso livello, per semplice lettura dei tipi piu' comuni */
   LowParser LowP;
     
   /* Stream in ingresso */
   InputStream* pIn;
   ifstream* pf;
   
   /* Buffer per le stringhe */
   char sStringBuf[iBufSize];
   char sStringBufWithSpaces[iBufSize];
      
   /* Parser delle espressioni matematiche, 
    * usato per acquisire valori sicuramente numerici */
   MathParser& MathP;
      
   /* Tabella dei simboli da riconoscere; puo' essere cambiata */
   KeyTable KeyT;   
   
   /* Token di basso ed alto livello */
   LowParser::Token CurrLowToken;
   Token CurrToken;
   
   virtual HighParser::Token FirstToken(void);
   virtual void NextToken(const char* sFuncName);
   
   flag fIsDescription(void);
   int iGetDescription_(const char* const s);
   void Set_(void);
   void Remark_(void);
   
 public:   
   HighParser(MathParser& MP, KeyTable& KT, InputStream& streamIn);
   virtual ~HighParser(void);

   virtual void PutKeyTable(KeyTable& KT);            /* Attacca una nuova KeyTable */
   virtual int GetLineNumber(void) const;             /* Numero di linea corrente */   
   virtual HighParser::ErrOut GetLineData(void) const;
   
   virtual MathParser& GetMathParser(void);           /* Restituisce il math parser */
   virtual void Close(void);                          /* "Chiude" i flussi */
   
   virtual int GetDescription(void);                  /* Legge una parola chiave */
   virtual void ExpectDescription(void);              /* si attende una descrizione */
   virtual void ExpectArg(void);                      /* si attende una lista di argomenti */

   virtual flag IsKeyWord(const char* sKeyWord);      /* 1 se trova la keyword sKeyWord */
   virtual int IsKeyWord(void);                       /* numero della keyword trovata */
   
   virtual flag fIsArg(void);                         /* 1 se e' atteso un argomento */
   virtual void PutBackSemicolon(void);               /* Se ha letto un ";" lo rimette a posto */
   virtual integer GetInt(int iDefval = 0);           /* legge un intero con il mathpar */
   virtual doublereal GetReal(double dDefval = 0.0);  /* legge un reale col mathpar */
   virtual int GetWord(void);                         /* legge una keyword */
   virtual const char* GetString(void);               /* legge una stringa */
   virtual const char* GetStringWithDelims(enum Delims Del = DEFAULTDELIM); 
                                                      /* stringa delimitata */

   virtual Vec3 GetVec3(void);                   /* vettore Vec3 */
   virtual Vec3 GetVec3(const Vec3& vDef);       /* vettore Vec3 */
   virtual Mat3x3 GetMatR2vec(void);             /* matrice R mediante i due vettori */
   virtual Mat3x3 GetMat3x3Sym(void);            /* matrice 3x3 simmetrica */
   virtual Mat3x3 GetMat3x3(void);               /* matrice 3x3 arbitraria */
   virtual Mat3x3 GetMat3x3(const Mat3x3& mDef);
   
   virtual Vec6 GetVec6(void);
   virtual Vec6 GetVec6(const Vec6& vDef);
   virtual Mat6x6 GetMat6x6(void);
   virtual Mat6x6 GetMat6x6(const Mat6x6& mDef);
   
   virtual inline doublereal Get(const doublereal& d) { return GetReal(double(d)); };
   virtual inline Vec3 Get(const Vec3& v) { return GetVec3(v); };
   virtual inline Mat3x3 Get(const Mat3x3& m) { return GetMat3x3(m); };
   virtual inline Vec6 Get(const Vec6& v) { return GetVec6(v); };
   virtual inline Mat6x6 Get(const Mat6x6& m) { return GetMat6x6(m); };
   
   virtual void GetMat6xN(Mat3xN& m1, Mat3xN& m2, integer iNumCols);
};

/* Le funzioni:
 *   ExpectDescription()
 *   ExpectArg()
 * informano il parser di cio' che e' atteso; di default il costruttore
 * setta ExpectDescription().
 * 
 * Le funzioni:
 *   GetDescription()
 *   IsKeyWord()
 *   GetWord()
 * restituiscono un intero che corrisponde alla posizione occupata nella
 * KeyTable dalle parole corrispondenti, oppure -1 se la parola non e'
 * trovata. Si noti che IsKeyWord(), in caso di esito negativo, ripristina 
 * l'istream. Tutte preparano poi il parser per la lettura successiva.
 * 
 * La funzione 
 *   IsKeyWord(const char*) 
 * restituisce 0 se non trova la parola e ripristina l'istream, altrimenti 
 * restituisce 1 e prepara il parser alla lettura successiva.
 * 
 * Le funzioni 
 *   GetInt(), 
 *   GetReal(), 
 *   GetString(), 
 *   GetStringWithDelims(enum Delims)
 * restituiscono i valori attesi e preparano il prser alla lettura successiva.
 */

/* HighParser - end */
   
extern ostream& operator << (ostream& out, const HighParser::ErrOut& err);

#endif
