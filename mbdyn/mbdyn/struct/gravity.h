/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2004
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

/* Accelerazione di gravita'
 * 
 * Elemento Gravity: contiene direzione e modulo, espresso mediante un driver,
 * dell'accelerazione di gravita'. E' un elemento unico (ne puo' essere
 * dichiarato uno solo) ed e' puntato da tutti gli elementi della classe
 * ElemGravityOwner, ovvero elementi che generano forze di inerzia
 * (per ora: Body, Beam).
 * 
 * Vi e' poi la classe GravityOwner, che contiene il puntatore all'elemento
 * Gravity. Da essa e' derivata la classe ElemGravityOwner. Quando l'elemento
 * viene costruito il puntatore e' nullo. Al termine della generazione
 * degli elementi, se e' definito l'elemento Gravity, tutti gli elementi
 * ElemGravityOwner vengono inizializzati con il puntatore all'elemento 
 * Gravity.
 * 
 * Si e' scelta la soluzione di un elemento per contenere questi dati
 * perche' in questo modo si acquista in generalita'. Infatti e' possibile
 * dare una dinamica all'accelerazione (in vista della generalizzazione del 
 * tipo di elemento) mediante l'aggiunta di gradi di liberta', ecc.
 * 
 * L'accelerazione e' ottenuta mediante la chiamata della funzione propria 
 * flag fGetAcceleration(Vec3&) da parte degli elementi ElemGravityOwner.
 * Il flag dice se e' definita l'accelerazione. 
 * In caso positivo, viene copiata nel vettore passato per reference.
 */

#ifndef GRAVITY_H
#define GRAVITY_H

#include <elem.h>
#include <tpldrive.h>

/* Gravity - begin */

class Gravity : public Elem, public TplDriveOwner<Vec3> {
 protected:
   Vec3 Acc;
   
 public:
   Gravity(const TplDriveCaller<Vec3>* pDC, flag fOut);
   
   virtual ~Gravity(void);
   
   virtual inline void* pGet(void) const { 
      return (void*)this;
   };
   
   /* Scrive il contributo dell'elemento al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const;
   
   /* Tipo dell'elemento (usato solo per debug ecc.) */
   virtual Elem::Type GetElemType(void) const { 
      return Elem::GRAVITY; 
   };
      
   /* funzioni di servizio */

   /* Il metodo iGetNumDof() serve a ritornare il numero di gradi di liberta'
    * propri che l'elemento definisce. Non e' virtuale in quanto serve a 
    * ritornare 0 per gli elementi che non possiedono gradi di liberta'.
    * Viene usato nella costruzione dei DofOwner e quindi deve essere 
    * indipendente da essi. In genere non comporta overhead in quanto il 
    * numero di dof aggiunti da un tipo e' una costante e non richede dati 
    * propri.
    * Il metodo pGetDofOwner() ritorna il puntatore al DofOwner dell'oggetto.
    * E' usato da tutti quelli che agiscono direttamente sui DofOwner.
    * Non e' virtuale in quanto ritorna NULL per tutti i tipi che non hanno
    * dof propri.
    * Il metodo GetDofType() ritorna, per ogni dof dell'elemento, l'ordine.
    * E' usato per completare i singoli Dof relativi all'elemento.
    */
   
   /* funzioni proprie */
   
   /* Dimensioni del workspace */
   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
      *piNumRows = 0;
      *piNumCols = 0;
   };
   
   /* assemblaggio jacobiano */
   virtual VariableSubMatrixHandler& 
     AssJac(VariableSubMatrixHandler& WorkMat,
	    doublereal dCoef, 
	    const VectorHandler& XCurr,
	    const VectorHandler& XPrimeCurr);
   
   /* assemblaggio residuo */
   virtual SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
				    doublereal dCoef,
				    const VectorHandler& XCurr, 
				    const VectorHandler& XPrimeCurr);

   virtual const Vec3& GetAcceleration(const Vec3& /* X */ ) const { 
      return Acc;
   };
};

/* Gravity - end */


/* GravityOwner - begin */

/* Classe base di elementi che generano forze di inerzia */

class GravityOwner {
 protected:
   const Gravity* pGravity;
   
 public:
   GravityOwner(void);
   virtual ~GravityOwner(void);

   void PutGravity(const Gravity* pG);
   virtual flag fGetAcceleration(const Vec3& X, Vec3& Acc) const;
};

/* GravityOwner - end */


/* ElemGravityOwner - begin */

class ElemGravityOwner : virtual public Elem, public GravityOwner {
 protected:
   /*
    * momento statico e momento di inerzia nel sistema globale
    */
   virtual Vec3 GetS_int(void) const {
      std::cerr << "warning: using default GetS_int()" << std::endl;
      return Vec3(0.);
   };

   virtual Mat3x3 GetJ_int(void) const {
      std::cerr << "warning: using default GetJ_int()" << std::endl;
      return Mat3x3(0.);
   };

 public:
   ElemGravityOwner(unsigned int uL, Elem::Type T, flag fOut);
   ~ElemGravityOwner(void);
   
   /* Usata per inizializzare la quantita' di moto */
   virtual void SetValue(VectorHandler& X, VectorHandler& XP) const = 0;

   /* Consente di effettuare un casting sicuro da Elem* a ElemGravityOwner* */
   virtual ElemGravityOwner* pGetElemGravityOwner(void) const { 
      return (ElemGravityOwner*)this; 
   };

   /*
    * massa
    */
   virtual doublereal dGetM(void) const {
      return 0.;
   };

   /*
    * momento statico e momento di inerzia trasportati nel punto X
    * e ruotati di R
    */
#if 0
   Vec3 GetS(const Vec3& X, const Mat3x3& R) const {
      return R*(GetS_int()-X*dGetM());
   };

   Mat3x3 GetJ(const Vec3& X, const Mat3x3& R) const {
      Vec3 S = X*dGetM()+GetS_int();
      return R*(GetJ_int()-Mat3x3(S, S/dGetM()))*R.Transpose();
   };
#endif
   
   Vec3 GetS(void) const {
	   return GetS_int();
   };
   
   Mat3x3 GetJ(void) const {
      return GetJ_int();
   };
   
#ifdef DEBUG
   virtual flag fIsElemGravityOwner(void) const { 
      return flag(1);
   };
#endif   
};

/* ElemGravityOwner - end */

#endif
