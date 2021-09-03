/*----------------------------------------------------------------------------*/
/* Copyright: CEA
 * contributors: F. Ledoux and N. Le Goff (2015)
 *
 * franck.ledoux@cea.fr
 * nicolas.le-goff@cea.fr
 *
 * The GMDS library is a computer program whose purpose is to provide a set of
 * functionnalities to represent and handle any type of meshes (2D, 3D,
 * triangles, tetrahedra, quad, hexa, polygons, polyhedra, etc.) and write
 * meshing algorithms. So it gathers many mathematical objects like points,
 * segment, quaternions, etc. and basic algorithms useful to build more evolved
 * ones.
 *
 * This software is governed by the CeCILL-C license under French law and
 * abiding by the rules of distribution of free software.  You can  use,
 * modify and/ or redistribute the software under the terms of the CeCILL-C
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 *
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability.
 *
 * In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or
 * data to be ensured and,  more generally, to use and operate it in the
 * same conditions as regards security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 */
/*
 * VariableManager.h
 *
 *  Created on: 03/10/2017
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef KMDS_VARIABLEMANAGER_H_
#define KMDS_VARIABLEMANAGER_H_
/*----------------------------------------------------------------------------*/
#include <string>
/*----------------------------------------------------------------------------*/
#include <KM/Utils/Exception.h>
#include <KM/Utils/KTypes.h>
#include <KM/Utils/Variable.h>
/*----------------------------------------------------------------------------*/
namespace kmds {
/*----------------------------------------------------------------------------*/
/** \class VariableManager
 *  \brief Handle the creation and update of a collection of variables. A
 *  	   variable is defined as a set of discrete values associated to a key.
 *  	   Few holes are in the key numerotation.
 */
class EXPORT_KMDS VariableManager
{
 public:
        /*------------------------------------------------------------------------*/
        /** \brief  Constructor.
        */
        VariableManager();

        /*------------------------------------------------------------------------*/
        /** \brief  Destructor.
        */
        ~VariableManager();

        /*------------------------------------------------------------------------*/
        /** \brief  creation of a variable allocated in the stack. The domain of the
         * 			variable is initialized to [0,initSize].
         */
        template <typename T>
        Variable<T>* createVariable(const T ADefaultValue, const std::string& AName, const int initSize = 2);

        /*------------------------------------------------------------------------*/
        /** \brief  Returns whether a variable exists.
        */
        bool doesVariableExist(const std::string& AName);

        /*------------------------------------------------------------------------*/
        /** \brief  Access to a variable.
        */
        template <typename T>
        EXPORT_KMDS Variable<T>* getVariable(const std::string& AName);

        /*------------------------------------------------------------------------*/
        /** \brief  suppression of a variable. the memory used in the stack is free.
        */
        void deleteVariable(const std::string& AName);

        /*------------------------------------------------------------------------*/
        /** \brief  suppression of a variable. the memory used in the stack is free.
        */
        void deleteVariable(VariableItf* AVar);

        /*------------------------------------------------------------------------*/
        /** \brief  set the domain of all the variables to \p ASize.
         */
        void resize(const int ASize);

        /*------------------------------------------------------------------------*/
        /** \brief  get the domain size of all the variables
                 */
        int getSize() const;

        /*------------------------------------------------------------------------*/
        /** \brief  get access to the i^th variable in a abstract form
         */
        VariableItf* getVariable(const TInt32 i);

        /*------------------------------------------------------------------------*/
        /** \brief  Initialize all the variables related to index \p AI. Default
         *          value is known by every single variable
         */
        void initializeVariables(const TCellID AI);
        /*------------------------------------------------------------------------*/
        /** \brief  Initialize all the variables related to indices \p AI to
         *          \p AI + \p ANb. Default value is known by every single variable
         */
        void initializeVariables(const TCellID AI, const TInt32 ANb);

        /*------------------------------------------------------------------------*/
        /** \brief indicates if variables are attached
         */
        bool empty() const;

        /*------------------------------------------------------------------------*/
        /** \brief compact all the variables
         */
        void compact();

        /*------------------------------------------------------------------------*/
        /** \brief serialize (*this) in AStr
         *
         * \param AStr an output streammap
         */
        void serialize(std::ostream& AStr);

        /*------------------------------------------------------------------------*/
        /** \brief unserialize (*this) from AStr
         *
         * \param AStr an input stream
         */
        void unserialize(std::istream& AStr);

        /*------------------------------------------------------------------------*/
        /** \brief  get the the number of variables
        */
        int getNbVariables() const;
        /*------------------------------------------------------------------------*/
        /** \brief  get the list of variables
        */
        std::vector<VariableItf*> getAllVariables();

 private:
        std::vector<VariableItf*> m_variables;
};
/*----------------------------------------------------------------------------*/
template <typename T>
Variable<T>*
VariableManager::createVariable(const T ADefaultValue, const std::string& AName, const int initSize)
{
        for (auto v : m_variables)
                if (v->getName() == AName) {
                        std::string mess = "Impossible to create a variable " + AName + ": name already used";
                        throw KException(mess);
                }

        Variable<T>* v = new Variable<T>(ADefaultValue, AName);
        m_variables.push_back(v);

        v->resize(initSize + 1);

        return v;
}
/*----------------------------------------------------------------------------*/
template <typename T>
Variable<T>*
VariableManager::getVariable(const std::string& AName)
{
        for (auto v : m_variables) {
                if (v->getName() == AName)
                        return dynamic_cast<Variable<T>*>(v);
        }
        std::string mess = "No variable named " + AName;
        throw KException(mess);
}
/*----------------------------------------------------------------------------*/
}  // namespace kmds
/*----------------------------------------------------------------------------*/
#endif /* KMDS_VARIABLEMANAGER_H_ */
/*----------------------------------------------------------------------------*/
