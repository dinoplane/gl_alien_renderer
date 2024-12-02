#include <cloth_system.hpp>
#include <shader_c.hpp>
#include <shader_s.hpp>
#include <gl_bindings.h>
#include <chrono>
#include <iostream>
#include<Eigen/Geometry>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include<Eigen/SparseCholesky>
//#include<Eigen/SparseLU>
#include<Eigen/IterativeLinearSolvers>
#include <pardiso.h>

//#include "pardiso.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>

static void dumpCSR(const char* filename, int n, int* ia, int* ja, double* a)
{
  std::fstream fout(filename, std::ios::out);
  fout << n << std::endl;
  fout << n << std::endl;
  fout << ia[n] << std::endl;
  
  for (int i = 0; i <= n; i++)
  {
    fout << ia[i] << std::endl;
  }

  for (int i = 0; i < ia[n]; i++)
  {
    fout << ja[i] << std::endl;
  }

  for (int i = 0; i < ia[n]; i++)
  {
    fout << a[i] << std::endl;
  }

  fout.close();
}


static void printCSR(int n, int nnz, int* ia, int* ja, double* a)
{
  std::cout << "rows: " << std::setw(10) << n << std::endl;
  std::cout << "nnz : " << std::setw(10) << nnz << std::endl;

  if (nnz == n*n)
  {
    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n; j++)
      {
        std::cout << std::setw(10) << a[i*n + j];
      }
      std::cout << std::endl;
    }
  }
  else
  {
    for (int i = 0; i < n; i++)
    {
      for (int index = ia[i]; index < ia[i+1]; index++)
      {
        int j = ja[index];
        std::cout << std::setw(10) << "(" << i << ", " << j << ") " << a[index];
      }
      std::cout << std::endl;
    }
  }
}


static void shiftIndices(int n, int nonzeros, int* ia, int* ja, int value)
{
  int i;
    for (i = 0; i < n+1; i++) 
    {
        ia[i] += value;
    }
    for (i = 0; i < nonzeros; i++) 
    {
        ja[i] += value;
    }
}




#define EIGEN_DONT_PARALLELIZE

static void doPardiso() {
/* Matrix data.
     */

    int    n =10;
    int    ia[11] = { 0, 6, 11, 15, 19, 22, 26, 29, 32, 33, 34};
    int    ja[34] = { 0, 2, 5, 6, 8, 9, 1, 2, 4, 8, 9, 2, 7, 8, 
                      9, 3, 6, 8, 9, 4, 8, 9, 5, 7, 8, 9, 6, 8, 
                      9, 7, 8, 9, 8, 9};
    
    double  a[34] = { 7.0,      1.0,           2.0, 7.0,     17.0, 8.5,
                          -4.0, 8.0,      2.0,                6.0, 3.0,
                                1.0,                     5.0, 6.0, 3.0,
                                     7.0,           9.0,     16.0, 8.0,
                                          0.0,               -4.0,-2.0,
                                               3.0,      8.0,18.0, 9.0,
                                                   11.0,     12.0, 6.0,
                                                         5.0, 4.0, 2.0,
                                                              0.0,
                                                                   0.0};



    int      nnz = ia[n];
    int      mtype = -2;        /* Real unsymmetric matrix */

    dumpCSR("matrix_sample_mtype_2.csr", n, ia, ja, a);
    printCSR(n, nnz, ia, ja, a);

    /* Internal solver memory pointer pt,                  */
    /* 32-bit: int pt[64]; 64-bit: long int pt[64]         */
    /* or void *pt[64] should be OK on both architectures  */ 
    void    *pt[64];

    /* Pardiso control parameters. */
    int      iparm[65];
    double   dparm[64];
    int      solver;
    int      maxfct, mnum, phase, error, msglvl;

    /* Number of processors. */
    int      num_procs;

    /* Auxiliary variables. */
    char    *var;
    int      i;

    double   ddum;              /* Double dummy */
    int      idum;              /* Integer dummy. */

/* -------------------------------------------------------------------- */
/* ..  Setup Pardiso control parameters and initialize the solvers      */
/*     internal adress pointers. This is only necessary for the FIRST   */
/*     call of the PARDISO solver.                                      */
/* ---------------------------------------------------------------------*/
      
    error = 0;
    solver = 0; /* use sparse direct solver */
    pardisoinit_d(pt,  &mtype, &solver, &iparm[1], dparm, &error);

    if (error != 0)
    {
        if (error == -10 )
           printf("No license file found \n");
        if (error == -11 )
           printf("License is expired \n");
        if (error == -12 )
           printf("Wrong username or hostname \n");
         
    }
    else
        printf("[PARDISO]: License check was successful ... \n");
 

    /* Numbers of processors, value of OMP_NUM_THREADS */
    // var = getenv("OMP_NUM_THREADS");
    // if(var != NULL)
    //     sscanf( var, "%d", &num_procs );
    // else 
    // {
    //     printf("Set environment OMP_NUM_THREADS to 1");
    //     // exit(1);
    // }
    num_procs = 12;
    iparm[3]  = num_procs;
    iparm[11] = 1;
    iparm[13] = 0;
   
    
    maxfct = 1;         /* Maximum number of numerical factorizations.  */
    mnum   = 1;         /* Which factorization to use. */
    
    msglvl = 1;         /* Print statistical information  */
    error  = 0;         /* Initialize error flag */


/* -------------------------------------------------------------------- */    
/* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
/*     notation.                                                        */
/* -------------------------------------------------------------------- */ 
    shiftIndices(n, nnz, ia, ja, 1);


/* -------------------------------------------------------------------- */
/*  .. pardiso_chk_matrix(...)                                          */
/*     Checks the consistency of the given matrix.                      */
/*     Use this functionality only for debugging purposes               */
/* -------------------------------------------------------------------- */
     pardiso_chkmatrix_d(&mtype, &n, a, ia, ja, &error);
     if (error != 0) 
     {
         printf("\nERROR in consistency of matrix: %d", error);
         // exit(1);
     }

 
// /* -------------------------------------------------------------------- */    
// /* ..  Reordering and Symbolic Factorization.  This step also allocates */
// /*     all memory that is necessary for the factorization.              */
// /* -------------------------------------------------------------------- */ 
    
    int nrows_S = 2;
    phase       = 12;
    iparm[38]   = nrows_S; 

    int nb = 0;

    pardiso_d(pt, &maxfct, &mnum, &mtype, &phase,
               &n, a, ia, ja, &idum, &nb,
               &iparm[1], &msglvl, &ddum, &ddum, &error,  dparm);
    
    if (error != 0) 
    {
        printf("\nERROR during symbolic factorization: %d", error);
        return;
    }
    printf("\nReordering completed ...\n");
    printf("Number of nonzeros in factors  = %d\n", iparm[18]);
    printf("Number of factorization GFLOPS = %d\n", iparm[19]);
    printf("Number of nonzeros is   S      = %d\n", iparm[39]);
   
/* -------------------------------------------------------------------- */    
/* ..  allocate memory for the Schur-complement and copy it there.      */
/* -------------------------------------------------------------------- */ 
    int nonzeros_S = iparm[39];
 
    int* iS        = new int[nrows_S+1];
    int* jS        = new int[nonzeros_S];
    double* S      = new double[nonzeros_S];
  
    pardiso_get_schur_d(pt, &maxfct, &mnum, &mtype, S, iS, jS);


/* -------------------------------------------------------------------- */    
/* ..  Convert matrix from 1-based Fortan notation to 0-based C         */
/*     notation.                                                        */
/* -------------------------------------------------------------------- */ 
    shiftIndices(n, nnz, ia, ja, -1);
    
/* -------------------------------------------------------------------- */    
/* ..  Convert Schur complement from Fortran notation held internally   */
/*     to 0-based C notation                                            */
/* -------------------------------------------------------------------- */ 
    shiftIndices(nrows_S, nonzeros_S, iS, jS, -1);

    printCSR(nrows_S, nonzeros_S, iS, jS, S);

/* -------------------------------------------------------------------- */    
/* ..  Termination and release of memory.                               */
/* -------------------------------------------------------------------- */ 
    phase = -1;                 /* Release internal memory. */

    pardiso_d(pt, &maxfct, &mnum, &mtype, &phase,
            &n, &ddum, ia, ja, &idum, &idum,
            &iparm[1], &msglvl, &ddum, &ddum, &error,  dparm);


    delete[] iS;
    delete[] jS;
    delete[] S;

    printf ("EXIT: Completed\n");

}


/*
def signedAngle(u = None,v = None,n = None):
    # This function calculates the signed angle between two vectors, "u" and "v",
    # using an optional axis vector "n" to determine the direction of the angle.
    #
    # Parameters:
    #   u: numpy array-like, shape (3,), the first vector.
    #   v: numpy array-like, shape (3,), the second vector.
    #   n: numpy array-like, shape (3,), the axis vector that defines the plane
    #      in which the angle is measured. It determines the sign of the angle.
    #
    # Returns:
    #   angle: double, the signed angle (in radians) from vector "u" to vector "v".
    #          The angle is positive if the rotation from "u" to "v" follows
    #          the right-hand rule with respect to the axis "n", and negative otherwise.
    #
    # The function works by:
    # 1. Computing the cross product "w" of "u" and "v" to find the vector orthogonal
    #    to both "u" and "v".
    # 2. Calculating the angle between "u" and "v" using the arctan2 function, which
    #    returns the angle based on the norm of "w" (magnitude of the cross product)
    #    and the dot product of "u" and "v".
    # 3. Using the dot product of "n" and "w" to determine the sign of the angle.
    #    If this dot product is negative, the angle is adjusted to be negative.
    #
    # Example:
    #   signedAngle(np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1]))
    #   This would return a positive angle (π/2 radians), as the rotation
    #   from the x-axis to the y-axis is counterclockwise when viewed along the z-axis.
    w = np.cross(u,v)
    angle = np.arctan2( np.linalg.norm(w), np.dot(u,v) )
    if (np.dot(n,w) < 0):
        angle = - angle

    return angle

def mmt(matrix):
    return matrix + matrix.T
*/

static double signedAngle(
    const Eigen::Ref<const Eigen::Vector3d>& u, 
    const Eigen::Ref<const Eigen::Vector3d>& v, 
    const Eigen::Ref<const Eigen::Vector3d>& n){
    Eigen::Vector3d w = u.cross(v);
    double angle = atan2(w.norm(), u.dot(v));
    if (n.dot(w) < 0){
        angle = -angle;
    }
    return angle;
}


static Eigen::Matrix3d mmt(
    const Eigen::Ref<const Eigen::Matrix3d>& matrix){
    return matrix + matrix.transpose();
}

/*
"""# Hinge angle, its gradient, and Hessian"""

#          x2
#          /\
#         /  \
#      e1/    \e3
#       /  t0  \
#      /        \
#     /    e0    \
#   x0------------x1
#     \          /
#      \   t1   /
#       \      /
#      e2\    /e4
#         \  /
#          \/
#          x3
#
#  Edge orientation: e0,e1,e2 point away from x0
#                       e3,e4 point away from x1

def getTheta(x0, x1 = None, x2 = None, x3 = None):

    if np.size(x0) == 12:  # Allow another type of input where x0 contains all the info
      x1 = x0[3:6]
      x2 = x0[6:9]
      x3 = x0[9:12]
      x0 = x0[0:3]

    m_e0 = x1 - x0
    m_e1 = x2 - x0
    m_e2 = x3 - x0

    n0 = np.cross(m_e0, m_e1)
    n1 = np.cross(m_e2, m_e0)

    # Calculate the signed angle using the provided function
    theta = signedAngle(n0, n1, m_e0)

    return theta
    */

static double getTheta(
        const Eigen::Ref< const Eigen::Vector3d > & node0, 
        const Eigen::Ref< const Eigen::Vector3d > & node1, 
        const Eigen::Ref< const Eigen::Vector3d > & node2, 
        const Eigen::Ref< const Eigen::Vector3d > & node3){
    Eigen::Vector3d m_e0 = node1 - node0;
    Eigen::Vector3d m_e1 = node2 - node0;
    Eigen::Vector3d m_e2 = node3 - node0;

    Eigen::Vector3d n0 = m_e0.cross(m_e1);
    Eigen::Vector3d n1 = m_e2.cross(m_e0);

    return signedAngle(n0, n1, m_e0);
}

void ClothSystem::InitializeSystemData(void* params) {
    ClothSystemParameters* clothParams = static_cast<ClothSystemParameters*>(params);
    clothParams->particleCount = (clothParams->clothSideLength + 1) * (clothParams->clothSideLength + 1);
    particleSystemDataBlock.particleCount = clothParams->particleCount;
    particleCount = clothParams->particleCount;
    particleSystemDataBlock.timeStep = clothParams->timeStep;
    dofCount = 3 * particleCount;
    vecValCount = particleCount * 4;
    
    indiceCount = 6 * clothParams->clothSideLength * clothParams->clothSideLength;

    double Y = clothParams->youngModulus;
    double h = clothParams->thickness;
    
    particleSystemDataBlock.bendingStiffness = 2.0 / sqrt(3.0) * Y * pow(h, 3.0) / 12.0;
    particleSystemDataBlock.timeStep = clothParams->timeStep;
    particleSystemDataBlock.tolerance = particleSystemDataBlock.bendingStiffness / (0.01) * 1.0e-3;
    particleSystemDataBlock.gravityAccel = clothParams->gravityAccel;


    bendingStiffness = particleSystemDataBlock.bendingStiffness;    
    dt = particleSystemDataBlock.timeStep;
    tol = particleSystemDataBlock.tolerance;
    gravityAccel = particleSystemDataBlock.gravityAccel;
    totalTime = 0.0;
    type = "base";
};



void ClothSystem::InitializeBufferData(void* params) {
    ClothSystemParameters* clothParams = static_cast<ClothSystemParameters*>(params);
    clothParams->particleCount = (clothParams->clothSideLength+ 1) * (clothParams->clothSideLength + 1);
    // count particles
    clothPositionVec.resize(vecValCount);
    // velocityVec.resize(vecValCount);
    // forceVec.resize(vecValCount);
    // normalsVec.resize(clothParams->particleCount);
    particleDataVec.resize(clothParams->particleCount);

    lastdofPositions.resize(dofCount);
    dofPositions.resize(dofCount);
    dofVelocities.resize(dofCount);

    //std::cout << "HAI" << lastdofPositions.size() << std::endl;

    // tangentsVec.resize(clothParams->particleCount);
    // bitangentsVec.resize(clothParams->particleCount);

    //dofVelocities.push_back(3.0);
    glm::dvec4 meshStartPosition = glm::dvec4(0.0, 0.0, 0.0, 1.0);
    double yIncrement = 0.5 * sqrt(3.0) * clothParams->cellSideLength;

    uint32_t clothSideParticleCount = clothParams->clothSideLength + 1;
    double deltaMass = clothParams->totalMass / clothParams->particleCount;
    massVector.resize(dofCount, deltaMass);
    // massMatrix.resize(dofCount*dofCount, 0.0);
    // for (uint32_t i = 0; i < dofCount; ++i) {
    //     massMatrix[i * dofCount + i] = deltaMass;
    // }


    // freeindices 0-9
    fixedNodes.push_back(0);
    fixedNodes.push_back(1);
    fixedNodes.push_back(2);

    //fixedNodes.push_back(41);
    //fixedNodes.push_back(42);
    //fixedNodes.push_back(43);
    //fixedNodes.push_back(67);

    //fixedNodes.push_back(190);

    //fixedNodes.push_back(250);

    //fixedNodes.push_back(330);


    //fixedNodes.push_back(8);

    uint32_t currFixedIdx = 0;
    for (uint32_t node = 0; node < particleCount; ++node) {
        if (currFixedIdx < fixedNodes.size() && node == fixedNodes[currFixedIdx]){
            ++currFixedIdx;
            fixedIndices.push_back(3 * node);
            fixedIndices.push_back(3 * node + 1);
            fixedIndices.push_back(3 * node + 2);

            oldToNewIndiceMapping.push_back(UINT32_MAX);
            oldToNewIndiceMapping.push_back(UINT32_MAX);
            oldToNewIndiceMapping.push_back(UINT32_MAX);

        
        }
        else {
            freeIndices.push_back(3 * node);
            freeIndices.push_back(3 * node + 1);
            freeIndices.push_back(3 * node + 2);

            oldToNewIndiceMapping.push_back(3 * node     - currFixedIdx * 3);
            oldToNewIndiceMapping.push_back(3 * node + 1 - currFixedIdx * 3);
            oldToNewIndiceMapping.push_back(3 * node + 2 - currFixedIdx * 3);
        }
    }
    freeDOFCount = freeIndices.size();
    fixedNodesCount = fixedNodes.size();

    massMatrix = Eigen::SparseMatrix<double>(freeDOFCount, freeDOFCount);
    SparseEntries massEntries;
    massEntries.reserve(freeDOFCount);
    for (uint32_t i = 0; i < freeDOFCount; ++i) {
        massEntries.push_back(Eigen::Triplet<double>(i, i, deltaMass));
    }
    massMatrix.setFromTriplets(massEntries.begin(), massEntries.end());


    // for (uint32_t i = 0; i < clothParams->clothSideLength; ++i) {

    // }

    for (uint32_t rowIdx = 0; rowIdx < clothSideParticleCount; ++rowIdx) {
        glm::dvec4 rowStartPosition =
            meshStartPosition;

        rowStartPosition.z -= yIncrement * rowIdx;
        //rowStartPosition.y = 0.5*sin(rowIdx * clothParams->cellSideLength * 2.0);
        if (rowIdx % 2 == 1) {
            rowStartPosition.x += 0.5 * clothParams->cellSideLength;
        }

        glm::dvec4 currentPosition = rowStartPosition;
        for (uint32_t colIdx = 0; colIdx < clothSideParticleCount; ++colIdx) {
            //currentPosition.y = rowStartPosition.y + 0.5 * sin(colIdx * clothParams->cellSideLength * 2.0);
            uint32_t particleIdx = rowIdx * clothSideParticleCount + colIdx;
            for (uint32_t i = 0; i < 4; ++i) {
                clothPositionVec[particleIdx * 4 + i] = currentPosition[i];
                //clothVelocityVec[particleIdx * 4 + i] = 0.0;
                //clothForceVec[particleIdx * 4 + i] = 0.0;
            }
            for (uint32_t i = 0; i < 3; ++i) {
                lastdofPositions[particleIdx * 3 + i] = currentPosition[i];
                dofPositions[particleIdx * 3 + i] = currentPosition[i];
                dofVelocities[particleIdx * 3 + i] = 0.0;
            }

            // clothPositionVec[particleIdx] = currentPosition;
            // velocityVec[particleIdx] = glm::dvec4(0.0);
            // forceVec[particleIdx] = glm::dvec4(0.0);
            particleDataVec[particleIdx].mass = deltaMass;

            // normalsVec[particleIdx] = glm::dvec4(0.0);
            currentPosition.x += clothParams->cellSideLength;
        }
    }

    // Fill in the EBO, edges, and hinges
    uint32_t c = clothParams->clothSideLength;
    edgeCount = c * (3 * c + 2);
    hingeCount = c * (3 * c - 2);
    edgeVec.resize(edgeCount);
    hingeVec.resize(hingeCount);
    uint32_t edgeIdx = 0;
    uint32_t hingeIdx = 0;
    
    indicesVec.resize(6 * clothParams->clothSideLength * clothParams->clothSideLength);
    for (uint32_t rowIdx = 0; rowIdx < clothParams->clothSideLength; ++rowIdx) {
        if (rowIdx % 2 == 1) {
            edgeVec[edgeIdx] = glm::ivec2(rowIdx * clothSideParticleCount, (rowIdx + 1) * clothSideParticleCount);
            edgeIdx += 1;
        }   
        for (uint32_t colIdx = 0; colIdx < clothParams->clothSideLength; ++colIdx) {
            // ADC / ABD
            //    G---------C---------D---------E   
            //     \   ____/ \   ____/ \   ____/ \  
            //      \ /       \ /       \ /       \ 
            //       H---------A---------B---------F
            //
            //    ABC / BDC
            //       G---------C---------D---------E
            //      / \____   / \____   / \____   / 
            //     /       \ /       \ /       \ /  
            //    H---------A---------B---------F   

            //                          Edges = Horizontals + Verticals
            //    0---1---2             Edges = c * (c + 1) + (2 * c + 1) * c 
            //    |\  |\  |             c * [(c + 1) + (2 * c + 1)] = c * (3 * c + 2) 
            //    |  \|  \|             2 * 3 + 5 * 2 = 6 + 10 = 16
            //    3---4---5             Hinges = Diagonals + Verticals + Horizontals
            //    |  /|  /|               
            //    |/  |/  |             Hinges = c * c + (c - 1) * c + (c - 1) * c
            //    6---7---8             Hinges = c * ( 3 * c - 2 )
            //         
            

            uint32_t cIdx = rowIdx * clothSideParticleCount + colIdx;
            uint32_t dIdx = cIdx + 1;
            uint32_t aIdx = (rowIdx + 1) * clothSideParticleCount + colIdx;
            uint32_t bIdx = aIdx + 1;

            uint32_t indiceIdx = 6 * (rowIdx * clothParams->clothSideLength + colIdx);
            //fmt::print("\ni {}\n", indiceIdx);

            uint32_t topHingeIdx = cIdx - clothSideParticleCount;
            uint32_t bottomHingeIdx = aIdx;
            uint32_t rightHingeIdx = cIdx + 1;
            uint32_t leftHingeIdx = aIdx - 1;

            if (rowIdx % 2 == 0) {
                indicesVec[indiceIdx + 0] = aIdx;
                indicesVec[indiceIdx + 1] = dIdx;
                indicesVec[indiceIdx + 2] = cIdx;

                indicesVec[indiceIdx + 3] = aIdx;
                indicesVec[indiceIdx + 4] = bIdx;
                indicesVec[indiceIdx + 5] = dIdx;

                hingeVec[hingeIdx] = glm::ivec4(aIdx, dIdx, cIdx, bIdx);

                edgeVec[edgeIdx + 0] = glm::ivec2(aIdx, dIdx);
                edgeVec[edgeIdx + 1] = glm::ivec2(dIdx, cIdx);
                edgeVec[edgeIdx + 2] = glm::ivec2(cIdx, aIdx);
                
            }
            else {
                indicesVec[indiceIdx + 0] = aIdx;
                indicesVec[indiceIdx + 1] = bIdx;
                indicesVec[indiceIdx + 2] = cIdx;

                indicesVec[indiceIdx + 3] = bIdx;
                indicesVec[indiceIdx + 4] = dIdx;
                indicesVec[indiceIdx + 5] = cIdx;
                
                hingeVec[hingeIdx] = glm::ivec4(cIdx, bIdx, dIdx, aIdx);
            
                edgeVec[edgeIdx + 0] = glm::ivec2(bIdx, dIdx);
                edgeVec[edgeIdx + 1] = glm::ivec2(dIdx, cIdx);
                edgeVec[edgeIdx + 2] = glm::ivec2(cIdx, bIdx);

                topHingeIdx = dIdx - clothSideParticleCount;
                bottomHingeIdx = bIdx;
                rightHingeIdx = aIdx + 1;
                leftHingeIdx = cIdx - 1;
            }
            hingeIdx += 1;
            edgeIdx += 3;

            if (colIdx > 0){
                hingeVec[hingeIdx] = glm::ivec4(cIdx, aIdx, rightHingeIdx, leftHingeIdx);
                hingeIdx += 1;
            } 

            if (rowIdx > 0){
                hingeVec[hingeIdx] = glm::ivec4(cIdx, dIdx, topHingeIdx, bottomHingeIdx);
                hingeIdx += 1;
            }

        }

        uint32_t lastRowParticleIdx = rowIdx * clothSideParticleCount + clothParams->clothSideLength;
        if (rowIdx % 2 == 0) {
            edgeVec[edgeIdx] = glm::ivec2(lastRowParticleIdx, lastRowParticleIdx + clothSideParticleCount);
            edgeIdx += 1;
        }   


    }
    uint32_t cIdx = clothParams->clothSideLength * clothSideParticleCount;
    for (uint32_t colIdx = 0; colIdx < clothParams->clothSideLength; ++colIdx){
        edgeVec[edgeIdx] = glm::ivec2(cIdx, cIdx + 1);
        edgeIdx += 1;
        cIdx += 1;
    }

    undeformedEdgeLengthVec.resize(edgeCount);
    elasticStretchingVec.resize(edgeCount);
    
    const double Y = clothParams->youngModulus;
    const double h = clothParams->thickness;
    std::cout << std::endl;
    for (uint32_t edgeIdx = 0; edgeIdx < edgeCount; ++edgeIdx) {
        glm::ivec2 edge = edgeVec[edgeIdx];
        glm::dvec4 edgeStart (clothPositionVec[4 * edge[0]], clothPositionVec[4 * edge[0] + 1], clothPositionVec[4 * edge[0] + 2], 1.0);
        glm::dvec4 edgeEnd (clothPositionVec[4 * edge[1]], clothPositionVec[4 * edge[1] + 1], clothPositionVec[4 * edge[1] + 2], 1.0);
        glm::dvec4 edgeVec = edgeEnd - edgeStart;
        undeformedEdgeLengthVec[edgeIdx] = glm::length(edgeVec);
        elasticStretchingVec[edgeIdx] = 0.5 * sqrt(3.0) * Y * h * pow(undeformedEdgeLengthVec[edgeIdx], 2.0);
        //fmt::print("({}, {}), ",edge[0], edge[1]);
        
    }

    Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> bufferStrides(4, 1);
    Eigen::Map<
        Eigen::Matrix<
                double,
                Eigen::Dynamic, 1>
    > q(dofPositions.data(), dofCount, 1);
    std::cout << std::endl;
    thetaBarVec.resize(hingeCount);
    for (uint32_t hingeIdx = 0; hingeIdx < hingeCount; hingeIdx++){
        glm::ivec4 hinge = hingeVec[hingeIdx];
        Eigen::Vector3d x0 = q.segment<3>(3 * hinge[0]);  
        Eigen::Vector3d x1 = q.segment<3>(3 * hinge[1]);
        Eigen::Vector3d x2 = q.segment<3>(3 * hinge[2]);
        Eigen::Vector3d x3 = q.segment<3>(3 * hinge[3]);

        thetaBarVec[hingeIdx] = getTheta(x0, x1, x2, x3);
        //fmt::print("({}, {}, {}, {}), ", hinge[0], hinge[1], hinge[2], hinge[3]);

    }

    std::cout << std::endl;
    externalForcesVec.resize(dofCount);
    for (uint32_t particleIdx = 0; particleIdx < particleCount; ++particleIdx) {
        externalForcesVec[3 * particleIdx] = 0.0;
        externalForcesVec[3 * particleIdx + 1] = -gravityAccel * particleDataVec[particleIdx].mass;
        externalForcesVec[3 * particleIdx + 2] = 0.0;
        // glm::dvec4 normal = glm::dvec4(0.0, 1.0, 0.0, 0.0);
        // normalsVec[particleIdx] = normal;
        //fmt::print("{: .8f}, {: .8f}, {: .8f}, ", dofPositions[3 * particleIdx], dofPositions[3 * particleIdx + 1], dofPositions[3 * particleIdx + 2]);
    }

    std::cout << std::endl;
}

void ClothSystem::InitializeShaders(void* params) {
    ClothSystemParameters* clothParams = static_cast<ClothSystemParameters*>(params);
    std::string vertPath = "./resources/shader/" + clothParams->shaderName + "_particle.vert";
    std::string fragPath = "./resources/shader/" + clothParams->shaderName + "_particle.frag";
    std::string compPath = "./resources/shader/" + clothParams->shaderName + "_particle.comp";
    std::string edgeDebugVertPath = "./resources/shader/" + clothParams->shaderName + "_edge.vert";
    std::string edgeDebugFragPath = "./resources/shader/" + clothParams->shaderName + "_edge.frag";
    std::string hingeDebugVertPath = "./resources/shader/" + clothParams->shaderName + "_hinge.vert";
    std::string hingeDebugFragPath = "./resources/shader/" + clothParams->shaderName + "_hinge.frag";
    std::string fixedNodeVertPath = "./resources/shader/" + clothParams->shaderName + "_fixed.vert";
    std::string fixedNodeFragPath = "./resources/shader/" + clothParams->shaderName + "_fixed.frag";

    particleShader = new Shader(vertPath.c_str(), fragPath.c_str());
    particleComputeShader = new ComputeShader(compPath.c_str());

    edgeDebugShader = new Shader(edgeDebugVertPath.c_str(), edgeDebugFragPath.c_str());
    hingeDebugShader = new Shader(hingeDebugVertPath.c_str(), hingeDebugFragPath.c_str());
    fixedNodesDebugShader = new Shader(fixedNodeVertPath.c_str(), fixedNodeFragPath.c_str());
}

void ClothSystem::InitializeBuffers(){
    // BaseParticleSystem::InitializeBuffers();
    glCreateBuffers(1, &positionBuffer);
    glNamedBufferStorage(
        positionBuffer, 
        sizeof(float) * vecValCount, 
        clothPositionVec.data(),
         GL_DYNAMIC_STORAGE_BIT
    );

    // glCreateBuffers(1, &velocityBuffer);
    // glNamedBufferStorage(
    //     velocityBuffer, 
    //     sizeof(double) * vecValCount, 
    //     clothVelocityVec.data(),
    //      GL_DYNAMIC_STORAGE_BIT
    // );

    // glCreateBuffers(1, &forcesBuffer);
    // glNamedBufferStorage(
    //     forcesBuffer, 
    //     sizeof(double) * vecValCount, 
    //     clothForceVec.data(),
    //      GL_DYNAMIC_STORAGE_BIT
    // );
    

    glCreateBuffers(1, &EBO);
    glNamedBufferStorage(
        EBO, 
        sizeof(uint32_t) * indiceCount, 
        indicesVec.data(),
         GL_DYNAMIC_STORAGE_BIT
    );

    // glCreateBuffers(1, &particleDataBuffer);
    // glNamedBufferStorage(
    //     particleDataBuffer, 
    //     sizeof(ParticleDataBlock) * particleCount, 
    //     nullptr,
    //      GL_DYNAMIC_STORAGE_BIT
    // );

     glCreateBuffers(1, &particleSystemDataBuffer);
     glNamedBufferStorage(
         particleSystemDataBuffer, 
         sizeof(ClothSystemDataBlock), 
         &particleSystemDataBlock,
          GL_DYNAMIC_STORAGE_BIT
     );


    glCreateBuffers(1, &edgesBuffer);
    glNamedBufferStorage(
       edgesBuffer, 
       sizeof(glm::ivec2) * edgeCount, 
       edgeVec.data(),
        GL_DYNAMIC_STORAGE_BIT
    );

    glCreateBuffers(1, &hingesBuffer);
    glNamedBufferStorage(
       hingesBuffer, 
       sizeof(glm::ivec4) * hingeCount, 
       hingeVec.data(),
        GL_DYNAMIC_STORAGE_BIT
    );

    glCreateBuffers(1, &edgeDebugDataBuffer);
    glNamedBufferStorage(
       edgeDebugDataBuffer, 
       sizeof(glm::dvec4) * 4, 
       edgeDebugDataVec,
        GL_DYNAMIC_STORAGE_BIT
    );

    glCreateBuffers(1, &edgeDebugEBO);
    glNamedBufferStorage(
       edgeDebugEBO, 
       sizeof(uint32_t) * 6, 
       edgeDebugEBOData,
        GL_DYNAMIC_STORAGE_BIT
    );
    
    glCreateBuffers(1, &pickedDebugDataBuffer);
    glNamedBufferStorage(
       pickedDebugDataBuffer, 
       sizeof(PickedClothData), 
       &pickedClothData,
        GL_DYNAMIC_STORAGE_BIT
    );

    glCreateBuffers(1, &fixedNodesDataBuffer);
    glNamedBufferStorage(
       fixedNodesDataBuffer, 
       sizeof(uint32_t) * fixedNodes.size(), 
       fixedNodes.data(),
        GL_DYNAMIC_STORAGE_BIT
    );

    //https://stackoverflow.com/questions/12399422/how-to-set-linker-flags-for-openmp-in-cmakes-try-compile-function
    fmt::print("{}", fixedNodesCount);
    glCreateVertexArrays(1, &fixedNodesVAO);

    glEnableVertexArrayAttrib(fixedNodesVAO, POSITION_ATTRIB_LOC);
    
    // nodeIdx attribute
    glVertexArrayAttribIFormat(fixedNodesVAO, POSITION_ATTRIB_LOC, 1, GL_UNSIGNED_INT, 0);
    glVertexArrayAttribBinding(fixedNodesVAO, POSITION_ATTRIB_LOC, 0);
    // Eigen::setNbThreads(12);
    // doPardiso();
    InitializePardiso();

}


/*
def hessTheta(x0, x1 = None, x2 = None, x3 = None):

    if np.size(x0) == 12:  # Allow another type of input where x0 contains all the info
      x1 = x0[3:6]
      x2 = x0[6:9]
      x3 = x0[9:12]
      x0 = x0[0:3]

    m_e0 = x1 - x0
    m_e1 = x2 - x0
    m_e2 = x3 - x0
    m_e3 = x2 - x1
    m_e4 = x3 - x1

    m_cosA1 = np.dot(m_e0, m_e1) / (np.linalg.norm(m_e0) * np.linalg.norm(m_e1))
    m_cosA2 = np.dot(m_e0, m_e2) / (np.linalg.norm(m_e0) * np.linalg.norm(m_e2))
    m_cosA3 = -np.dot(m_e0, m_e3) / (np.linalg.norm(m_e0) * np.linalg.norm(m_e3))
    m_cosA4 = -np.dot(m_e0, m_e4) / (np.linalg.norm(m_e0) * np.linalg.norm(m_e4))

    m_sinA1 = np.linalg.norm(np.cross(m_e0, m_e1)) / (np.linalg.norm(m_e0) * np.linalg.norm(m_e1))
    m_sinA2 = np.linalg.norm(np.cross(m_e0, m_e2)) / (np.linalg.norm(m_e0) * np.linalg.norm(m_e2))
    m_sinA3 = -np.linalg.norm(np.cross(m_e0, m_e3)) / (np.linalg.norm(m_e0) * np.linalg.norm(m_e3))
    m_sinA4 = -np.linalg.norm(np.cross(m_e0, m_e4)) / (np.linalg.norm(m_e0) * np.linalg.norm(m_e4))

    m_nn1 = np.cross(m_e0, m_e3)
    m_nn1 /= np.linalg.norm(m_nn1)
    m_nn2 = -np.cross(m_e0, m_e4)
    m_nn2 /= np.linalg.norm(m_nn2)

    m_h1 = np.linalg.norm(m_e0) * m_sinA1
    m_h2 = np.linalg.norm(m_e0) * m_sinA2
    m_h3 = -np.linalg.norm(m_e0) * m_sinA3
    m_h4 = -np.linalg.norm(m_e0) * m_sinA4
    m_h01 = np.linalg.norm(m_e1) * m_sinA1
    m_h02 = np.linalg.norm(m_e2) * m_sinA2

    # Gradient of Theta (as an intermediate step)
    grad_theta = np.zeros((12, 1))
    grad_theta[0:3] = (m_cosA3 * m_nn1 / m_h3 + m_cosA4 * m_nn2 / m_h4).reshape(-1, 1)
    grad_theta[3:6] = (m_cosA1 * m_nn1 / m_h1 + m_cosA2 * m_nn2 / m_h2).reshape(-1, 1)
    grad_theta[6:9] = (-m_nn1 / m_h01).reshape(-1, 1)
    grad_theta[9:12] = (-m_nn2 / m_h02).reshape(-1, 1)

    # Intermediate matrices for Hessian
    m_m1 = np.cross(m_nn1, m_e1) / np.linalg.norm(m_e1)
    m_m2 = -np.cross(m_nn2, m_e2) / np.linalg.norm(m_e2)
    m_m3 = -np.cross(m_nn1, m_e3) / np.linalg.norm(m_e3)
    m_m4 = np.cross(m_nn2, m_e4) / np.linalg.norm(m_e4)
    m_m01 = -np.cross(m_nn1, m_e0) / np.linalg.norm(m_e0)
    m_m02 = np.cross(m_nn2, m_e0) / np.linalg.norm(m_e0)

    # Hessian matrix components
    M331 = m_cosA3 / (m_h3 ** 2) * np.outer(m_m3, m_nn1)
    M311 = m_cosA3 / (m_h3 * m_h1) * np.outer(m_m1, m_nn1)
    M131 = m_cosA1 / (m_h1 * m_h3) * np.outer(m_m3, m_nn1)
    M3011 = m_cosA3 / (m_h3 * m_h01) * np.outer(m_m01, m_nn1)
    M111 = m_cosA1 / (m_h1 ** 2) * np.outer(m_m1, m_nn1)
    M1011 = m_cosA1 / (m_h1 * m_h01) * np.outer(m_m01, m_nn1)

    M442 = m_cosA4 / (m_h4 ** 2) * np.outer(m_m4, m_nn2)
    M422 = m_cosA4 / (m_h4 * m_h2) * np.outer(m_m2, m_nn2)
    M242 = m_cosA2 / (m_h2 * m_h4) * np.outer(m_m4, m_nn2)
    M4022 = m_cosA4 / (m_h4 * m_h02) * np.outer(m_m02, m_nn2)
    M222 = m_cosA2 / (m_h2 ** 2) * np.outer(m_m2, m_nn2)
    M2022 = m_cosA2 / (m_h2 * m_h02) * np.outer(m_m02, m_nn2)

    B1 = 1 / np.linalg.norm(m_e0) ** 2 * np.outer(m_nn1, m_m01)
    B2 = 1 / np.linalg.norm(m_e0) ** 2 * np.outer(m_nn2, m_m02)

    N13 = 1 / (m_h01 * m_h3) * np.outer(m_nn1, m_m3)
    N24 = 1 / (m_h02 * m_h4) * np.outer(m_nn2, m_m4)
    N11 = 1 / (m_h01 * m_h1) * np.outer(m_nn1, m_m1)
    N22 = 1 / (m_h02 * m_h2) * np.outer(m_nn2, m_m2)
    N101 = 1 / (m_h01 ** 2) * np.outer(m_nn1, m_m01)
    N202 = 1 / (m_h02 ** 2) * np.outer(m_nn2, m_m02)

    # Initialize Hessian of Theta
    hess_theta = np.zeros((12, 12))

    hess_theta[0:3, 0:3] = mmt(M331) - B1 + mmt(M442) - B2
    hess_theta[0:3, 3:6] = M311 + M131.T + B1 + M422 + M242.T + B2
    hess_theta[0:3, 6:9] = M3011 - N13
    hess_theta[0:3, 9:12] = M4022 - N24
    hess_theta[3:6, 3:6] = mmt(M111) - B1 + mmt(M222) - B2
    hess_theta[3:6, 6:9] = M1011 - N11
    hess_theta[3:6, 9:12] = M2022 - N22
    hess_theta[6:9, 6:9] = -mmt(N101)
    hess_theta[9:12, 9:12] = -mmt(N202)

    # Make the Hessian symmetric
    hess_theta[3:6, 0:3] = hess_theta[0:3, 3:6].T
    hess_theta[6:9, 0:3] = hess_theta[0:3, 6:9].T
    hess_theta[9:12, 0:3] = hess_theta[0:3, 9:12].T
    hess_theta[6:9, 3:6] = hess_theta[3:6, 6:9].T
    hess_theta[9:12, 3:6] = hess_theta[3:6, 9:12].T

    return hess_theta
    */
/*
#          x2
#          /\
#         /  \
#      e1/    \e3
#       /  t0  \
#      /        \
#     /    e0    \
#   x0------------x1
#     \          /
#      \   t1   /
#       \      /
#      e2\    /e4
#         \  /
#          \/
#          x3
#
#  Edge orientation : e0, e1, e2 point away from x0
#                       e3, e4 point away from x1
*/

static Eigen::MatrixXd uvT(
    const Eigen::Ref< const Eigen::Vector3d> & u, 
    const Eigen::Ref< const Eigen::Vector3d> & v){
    return u * v.transpose();
}

static void hessTheta(
        const Eigen::Ref< const Eigen::Vector3d > & node0, 
        const Eigen::Ref< const Eigen::Vector3d > & node1, 
        const Eigen::Ref< const Eigen::Vector3d > & node2, 
        const Eigen::Ref< const Eigen::Vector3d > & node3, 
        double* hess){
    
    Eigen::Vector3d m_e0 = node1 - node0;
    Eigen::Vector3d m_e1 = node2 - node0;
    Eigen::Vector3d m_e2 = node3 - node0;
    Eigen::Vector3d m_e3 = node2 - node1;
    Eigen::Vector3d m_e4 = node3 - node1;

    double m_cosA1 = m_e0.dot(m_e1) / (m_e0.norm() * m_e1.norm());
    double m_cosA2 = m_e0.dot(m_e2) / (m_e0.norm() * m_e2.norm());
    double m_cosA3 = -m_e0.dot(m_e3) / (m_e0.norm() * m_e3.norm());
    double m_cosA4 = -m_e0.dot(m_e4) / (m_e0.norm() * m_e4.norm());

    double m_sinA1 = m_e0.cross(m_e1).norm() / (m_e0.norm() * m_e1.norm());
    double m_sinA2 = m_e0.cross(m_e2).norm() / (m_e0.norm() * m_e2.norm());
    double m_sinA3 = -m_e0.cross(m_e3).norm() / (m_e0.norm() * m_e3.norm());
    double m_sinA4 = -m_e0.cross(m_e4).norm() / (m_e0.norm() * m_e4.norm());

    Eigen::Vector3d m_nn1 = m_e0.cross(m_e3);
    m_nn1 = m_nn1 / m_nn1.norm();
    Eigen::Vector3d m_nn2 = -m_e0.cross(m_e4);
    m_nn2 = m_nn2 / m_nn2.norm();

    double m_h1 = m_e0.norm() * m_sinA1;
    double m_h2 = m_e0.norm() * m_sinA2;
    double m_h3 = -m_e0.norm() * m_sinA3;
    double m_h4 = -m_e0.norm() * m_sinA4;
    double m_h01 = m_e1.norm() * m_sinA1;
    double m_h02 = m_e2.norm() * m_sinA2;

    Eigen::Vector3d m_m1 = m_nn1.cross(m_e1) / m_e1.norm();
    Eigen::Vector3d m_m2 = -m_nn2.cross(m_e2) / m_e2.norm();
    Eigen::Vector3d m_m3 = -m_nn1.cross(m_e3) / m_e3.norm();
    Eigen::Vector3d m_m4 = m_nn2.cross(m_e4) / m_e4.norm();
    Eigen::Vector3d m_m01 = -m_nn1.cross(m_e0) / m_e0.norm();
    Eigen::Vector3d m_m02 = m_nn2.cross(m_e0) / m_e0.norm();
    //
    Eigen::Map<Eigen::MatrixXd> hess_theta(hess, 12, 12);
     hess_theta.setZero();

    Eigen::MatrixXd M331 = m_cosA3 / (m_h3 * m_h3) * uvT(m_m3, m_nn1);
    Eigen::MatrixXd M311 = m_cosA3 / (m_h3 * m_h1) * uvT(m_m1, m_nn1);
    Eigen::MatrixXd M131 = m_cosA1 / (m_h1 * m_h3) * uvT(m_m3, m_nn1);
    Eigen::MatrixXd M3011 = m_cosA3 / (m_h3 * m_h01) * uvT(m_m01, m_nn1);
    Eigen::MatrixXd M111 = m_cosA1 / (m_h1 * m_h1) * uvT(m_m1, m_nn1);
    Eigen::MatrixXd M1011 = m_cosA1 / (m_h1 * m_h01) * uvT(m_m01, m_nn1);
    
    Eigen::MatrixXd M442 = m_cosA4 / (m_h4 * m_h4) * uvT(m_m4, m_nn2);
    Eigen::MatrixXd M422 = m_cosA4 / (m_h4 * m_h2) * uvT(m_m2, m_nn2);
    Eigen::MatrixXd M242 = m_cosA2 / (m_h2 * m_h4) * uvT(m_m4, m_nn2);
    Eigen::MatrixXd M4022 = m_cosA4 / (m_h4 * m_h02) * uvT(m_m02, m_nn2);
    Eigen::MatrixXd M222 = m_cosA2 / (m_h2 * m_h2) * uvT(m_m2, m_nn2);
    Eigen::MatrixXd M2022 = m_cosA2 / (m_h2 * m_h02) * uvT(m_m02, m_nn2);

    Eigen::MatrixXd B1 = 1 / m_e0.norm() / m_e0.norm() * uvT(m_nn1, m_m01);
    Eigen::MatrixXd B2 = 1 / m_e0.norm() / m_e0.norm() * uvT(m_nn2, m_m02);

    Eigen::MatrixXd N13 = 1 / (m_h01 * m_h3) * uvT(m_nn1, m_m3);
    Eigen::MatrixXd N24 = 1 / (m_h02 * m_h4) * uvT(m_nn2, m_m4);
    Eigen::MatrixXd N11 = 1 / (m_h01 * m_h1) * uvT(m_nn1, m_m1);
    Eigen::MatrixXd N22 = 1 / (m_h02 * m_h2) * uvT(m_nn2, m_m2);
    Eigen::MatrixXd N101 = 1 / (m_h01 * m_h01) * uvT(m_nn1, m_m01);
    Eigen::MatrixXd N202 = 1 / (m_h02 * m_h02) * uvT(m_nn2, m_m02);

    //std::cout << M331 << std::endl;

    hess_theta.block<3, 3>(0, 0) = mmt(M331) - B1 + mmt(M442) - B2;
    hess_theta.block<3, 3>(0, 3) = M311 + M131.transpose() + B1 + M422 + M242.transpose() + B2;
    hess_theta.block<3, 3>(0, 6) = M3011 - N13;
    hess_theta.block<3, 3>(0, 9) = M4022 - N24;
    hess_theta.block<3, 3>(3, 3) = mmt(M111) - B1 + mmt(M222) - B2;
    hess_theta.block<3, 3>(3, 6) = M1011 - N11;
    hess_theta.block<3, 3>(3, 9) = M2022 - N22;
    hess_theta.block<3, 3>(6, 6) = -mmt(N101);
    hess_theta.block<3, 3>(9, 9) = -mmt(N202);

    hess_theta.block<3, 3>(3, 0) = hess_theta.block<3, 3>(0, 3).transpose();
    hess_theta.block<3, 3>(6, 0) = hess_theta.block<3, 3>(0, 6).transpose();
    hess_theta.block<3, 3>(9, 0) = hess_theta.block<3, 3>(0, 9).transpose();
    hess_theta.block<3, 3>(6, 3) = hess_theta.block<3, 3>(3, 6).transpose();
    hess_theta.block<3, 3>(9, 3) = hess_theta.block<3, 3>(3, 9).transpose();

}


    /*
# In the original code, there are probaly TWO sign errors in the expressions for m_h3 and m_h4.
# [Original code: % https://github.com/shift09/plates-shells/blob/master/src/bending.cpp]
# I indicated those two corrections by writing the word "CORRECTION" next
# to them.

def gradTheta(x0, x1 = None, x2 = None, x3 = None):

    if np.size(x0) == 12:  # Allow another type of input where x0 contains all the info
      x1 = x0[3:6]
      x2 = x0[6:9]
      x3 = x0[9:12]
      x0 = x0[0:3]

    m_e0 = x1 - x0
    m_e1 = x2 - x0
    m_e2 = x3 - x0
    m_e3 = x2 - x1
    m_e4 = x3 - x1

    m_cosA1 = np.dot(m_e0, m_e1) / (np.linalg.norm(m_e0) * np.linalg.norm(m_e1))
    m_cosA2 = np.dot(m_e0, m_e2) / (np.linalg.norm(m_e0) * np.linalg.norm(m_e2))
    m_cosA3 = -np.dot(m_e0, m_e3) / (np.linalg.norm(m_e0) * np.linalg.norm(m_e3))
    m_cosA4 = -np.dot(m_e0, m_e4) / (np.linalg.norm(m_e0) * np.linalg.norm(m_e4))

    m_sinA1 = np.linalg.norm(np.cross(m_e0, m_e1)) / (np.linalg.norm(m_e0) * np.linalg.norm(m_e1))
    m_sinA2 = np.linalg.norm(np.cross(m_e0, m_e2)) / (np.linalg.norm(m_e0) * np.linalg.norm(m_e2))
    m_sinA3 = -np.linalg.norm(np.cross(m_e0, m_e3)) / (np.linalg.norm(m_e0) * np.linalg.norm(m_e3))
    m_sinA4 = -np.linalg.norm(np.cross(m_e0, m_e4)) / (np.linalg.norm(m_e0) * np.linalg.norm(m_e4))

    m_nn1 = np.cross(m_e0, m_e3)
    m_nn1 = m_nn1 / np.linalg.norm(m_nn1)
    m_nn2 = -np.cross(m_e0, m_e4)
    m_nn2 = m_nn2 / np.linalg.norm(m_nn2)

    m_h1 = np.linalg.norm(m_e0) * m_sinA1
    m_h2 = np.linalg.norm(m_e0) * m_sinA2
    m_h3 = -np.linalg.norm(m_e0) * m_sinA3  # CORRECTION
    m_h4 = -np.linalg.norm(m_e0) * m_sinA4  # CORRECTION
    m_h01 = np.linalg.norm(m_e1) * m_sinA1
    m_h02 = np.linalg.norm(m_e2) * m_sinA2

    # Initialize the gradient
    gradTheta = np.zeros(12)

    gradTheta[0:3] = m_cosA3 * m_nn1 / m_h3 + m_cosA4 * m_nn2 / m_h4
    gradTheta[3:6] = m_cosA1 * m_nn1 / m_h1 + m_cosA2 * m_nn2 / m_h2
    gradTheta[6:9] = -m_nn1 / m_h01
    gradTheta[9:12] = -m_nn2 / m_h02

    return gradTheta
*/

static void gradTheta(
        const Eigen::Ref< const Eigen::Vector3d >& node0, 
        const Eigen::Ref< const Eigen::Vector3d >& node1, 
        const Eigen::Ref< const Eigen::Vector3d >& node2, 
        const Eigen::Ref< const Eigen::Vector3d >& node3, 
        double* grad){
    
    Eigen::Vector3d m_e0 = node1 - node0;
    Eigen::Vector3d m_e1 = node2 - node0;
    Eigen::Vector3d m_e2 = node3 - node0;
    Eigen::Vector3d m_e3 = node2 - node1;
    Eigen::Vector3d m_e4 = node3 - node1;

    double m_cosA1 = m_e0.dot(m_e1) / (m_e0.norm() * m_e1.norm());
    double m_cosA2 = m_e0.dot(m_e2) / (m_e0.norm() * m_e2.norm());
    double m_cosA3 = -m_e0.dot(m_e3) / (m_e0.norm() * m_e3.norm());
    double m_cosA4 = -m_e0.dot(m_e4) / (m_e0.norm() * m_e4.norm());

    double m_sinA1 = m_e0.cross(m_e1).norm() / (m_e0.norm() * m_e1.norm());
    double m_sinA2 = m_e0.cross(m_e2).norm() / (m_e0.norm() * m_e2.norm());
    double m_sinA3 = -m_e0.cross(m_e3).norm() / (m_e0.norm() * m_e3.norm());
    double m_sinA4 = -m_e0.cross(m_e4).norm() / (m_e0.norm() * m_e4.norm());

    Eigen::Vector3d m_nn1 = m_e0.cross(m_e3);
    m_nn1 /= m_nn1.norm();
    Eigen::Vector3d m_nn2 = -m_e0.cross(m_e4);
    m_nn2 /= m_nn2.norm();

    double m_h1 = m_e0.norm() * m_sinA1;
    double m_h2 = m_e0.norm() * m_sinA2;
    double m_h3 = -m_e0.norm() * m_sinA3;
    double m_h4 = -m_e0.norm() * m_sinA4;
    double m_h01 = m_e1.norm() * m_sinA1;
    double m_h02 = m_e2.norm() * m_sinA2;

    Eigen::Vector3d gradTheta0 = m_cosA3 * m_nn1 / m_h3 + m_cosA4 * m_nn2 / m_h4;
    Eigen::Vector3d gradTheta1 = m_cosA1 * m_nn1 / m_h1 + m_cosA2 * m_nn2 / m_h2;
    Eigen::Vector3d gradTheta2 = -m_nn1 / m_h01;
    Eigen::Vector3d gradTheta3 = -m_nn2 / m_h02;

    Eigen::Map<Eigen::Matrix<double, 12, 1>> dFMap(grad);
    dFMap.segment<3>(0) = gradTheta0;
    dFMap.segment<3>(3) = gradTheta1;
    dFMap.segment<3>(6) = gradTheta2;
    dFMap.segment<3>(9) = gradTheta3;


}

    /*
"""# Stretching energy for a shell, it's gradient, and Hessian"""

def gradEs_hessEs(node0 = None,node1 = None,l_k = None,EA = None):

# Inputs:
# node0: 1x3 vector - position of the first node
# node1: 1x3 vector - position of the last node

# l_k: reference length (undeformed) of the edge
# EA: scalar - stretching stiffness - Young's modulus times area

# Outputs:
# dF: 6x1  vector - gradient of the stretching energy between node0 and node 1.
# dJ: 6x6 vector - hessian of the stretching energy between node0 and node 1.

    ## Gradient of Es
    edge = node1 - node0

    edgeLen = np.linalg.norm(edge)
    tangent = edge / edgeLen
    epsX = edgeLen / l_k - 1
    dF_unit = EA * tangent * epsX
    dF = np.zeros((6))
    dF[0:3] = - dF_unit
    dF[3:6] = dF_unit

    ## Hessian of Es
    Id3 = np.eye(3)
    M = EA * ((1 / l_k - 1 / edgeLen) * Id3 + 1 / edgeLen * ( np.outer( edge, edge ) ) / edgeLen ** 2)

    dJ = np.zeros((6,6))
    dJ[0:3,0:3] = M
    dJ[3:6,3:6] = M
    dJ[0:3,3:6] = - M
    dJ[3:6,0:3] = - M
    return dF,dJ
*/



static void subFromSparseEntries(const Eigen::Ref< const Eigen::MatrixXd > & mat, 
                                uint32_t* indices, 
                                SparseEntries* sparseEntries,
                                const std::span<uint32_t> oldToNewIndices){
    for (uint32_t i = 0; i < mat.rows(); ++i){
        uint32_t oldRowIdx = indices[i];
        uint32_t newRowIdx = oldToNewIndices[oldRowIdx];
        if (newRowIdx == UINT32_MAX) {
            continue;
        }

        for (uint32_t j = 0; j < mat.cols(); ++j){
            uint32_t oldColIdx = indices[j];
            uint32_t newColIdx = oldToNewIndices[oldColIdx];
            if (newColIdx == UINT32_MAX) {
                continue;
            }

            if (mat(i, j) == 0.0) {
                continue;
            }
            
            sparseEntries->push_back(Eigen::Triplet<double>(newRowIdx, newColIdx, -mat(i, j)));
        }
    }
}

static void calculateGradEb_HessEb_Shell(
    const Eigen::Ref< const Eigen::Vector3d > & node0, 
    const Eigen::Ref< const Eigen::Vector3d > & node1, 
    const Eigen::Ref< const Eigen::Vector3d > & node2, 
    const Eigen::Ref< const Eigen::Vector3d > & node3, 
    double thetaBar, double kb, double* dF, double* dJ){
    /*
def gradEb_hessEb_Shell(x0, x1=None, x2=None, x3=None, theta_bar=0, kb=1.0):
    """
    Compute the gradient and Hessian of the bending energy for a shell.

    Parameters:
    x0 (array): Can either be a 3-element array (single point) or a 12-element array.
    x1, x2, x3 (arrays): Optional, 3-element arrays specifying points.
    theta_bar (double): Reference angle.
    kb (double): Bending stiffness.

    Returns:
    dF (array): Gradient of the bending energy.
    dJ (array): Hessian of the bending energy.
    """
    # Allow another type of input where x0 contains all the information
    if np.size(x0) == 12:
        x1 = x0[3:6]
        x2 = x0[6:9]
        x3 = x0[9:12]
        x0 = x0[:3]

    # Compute theta, gradient, and Hessian
    theta = getTheta(x0, x1, x2, x3)  # Replace with your getTheta function in Python
    grad = gradTheta(x0, x1, x2, x3)  # Replace with your gradTheta function in Python

    dF = 0.5 * kb * (2 * (theta - theta_bar) * grad)

    hess = hessTheta(x0, x1, x2, x3)  # Replace with your hessTheta function in Python
    dJ = 0.5 * kb * (2 * np.outer(grad, grad) + 2 * (theta - theta_bar) * hess)

    return dF, dJ

    */
    
    double theta = getTheta(node0, node1, node2, node3);
    double grad[12];
    gradTheta(node0, node1, node2, node3, grad);
    
    Eigen::Map<
        Eigen::Matrix<
            double,
            12, 1>
    > gradMap(grad);
    Eigen::Map<
        Eigen::Matrix<
            double,
            12, 1>
    > dFMap(dF);

    dFMap.segment<12>(0) = 0.5 * kb * (2.0 * (theta - thetaBar) * gradMap);
    
    double hess[144];
    hessTheta(node0, node1, node2, node3, hess);
    
    
    Eigen::Map<Eigen::Matrix<double, 12, 12>> dJMap(dJ);
    Eigen::Map<Eigen::Matrix<double, 12, 12>> hessMap(hess);
    //// TODO Which is the correct transpose?
    
    Eigen::MatrixXd gradgradT = gradMap*gradMap.transpose() ;
    dJMap.block<12, 12>(0, 0) = 0.5 * kb * (2.0 * gradgradT + 2.0 * (theta - thetaBar) * hessMap);
    
    

}

static void calculateFb_Jb(double* bendingForcesVecPtr, SparseEntries* sparseEntries, // double* bendingHessiansVecPtr, 
                        const std::vector<glm::ivec4>& hingeVec, 
                        const Eigen::Ref<const Eigen::Matrix<double, Eigen::Dynamic, 1>>& q, 
                        double kb, 
                        const std::vector<double>& thetaBar, 
                        uint32_t dofCount, uint32_t hingeCount,
                        const std::span<uint32_t> oldToNewIndiceMapping) {
    
    Eigen::Map<Eigen::Matrix<
            double,
            Eigen::Dynamic, 1
            >> bendingForces(bendingForcesVecPtr, dofCount, 1);
    // Eigen::Map<
    //     Eigen::Matrix<
    //     double,
    //     Eigen::Dynamic, Eigen::Dynamic
    //     >> bendingHessians(bendingHessiansVecPtr, dofCount, dofCount);


    // Bending Forces
    for (uint32_t hingeIdx = 0; hingeIdx < hingeCount; ++hingeIdx) {
        glm::ivec4 kHinge = hingeVec[hingeIdx];
        uint32_t node0 = kHinge[0];
        uint32_t node1 = kHinge[1];
        uint32_t node2 = kHinge[2];
        uint32_t node3 = kHinge[3];

        Eigen::Vector3d x0 = q.segment<3>(3 * node0);
        Eigen::Vector3d x1 = q.segment<3>(3 * node1);
        Eigen::Vector3d x2 = q.segment<3>(3 * node2);
        Eigen::Vector3d x3 = q.segment<3>(3 * node3);

        uint32_t ind[] = { 3 * node0, 3 * node0 + 1, 3 * node0 + 2,
                      3 * node1, 3 * node1 + 1, 3 * node1 + 2,
                      3 * node2, 3 * node2 + 1, 3 * node2 + 2,
                      3 * node3, 3 * node3 + 1, 3 * node3 + 2 };

        double dF[12];
        double dJ[144];


        Eigen::Map<
        Eigen::Matrix<
            double,
            12, 1
        >> dFMap(dF);

        Eigen::Map<
            Eigen::Matrix<
            double,
            12, 12>
        > dJMap(dJ);
        calculateGradEb_HessEb_Shell(x0, x1, x2, x3, thetaBar[hingeIdx], kb, dF, dJ);

        bendingForces(ind) = bendingForces(ind).eval() - dFMap;
        // bendingHessians(ind, ind) = bendingHessians(ind, ind).eval() - dJMap;
        subFromSparseEntries(dJMap, ind, sparseEntries, oldToNewIndiceMapping);
    }
}


static void calculateGradEs_HessEs(
    const Eigen::Ref< const Eigen::Vector3d>& node0, 
    const Eigen::Ref< const Eigen::Vector3d >& node1,
    double lk, double EA, double* dF, double* dJ) {
    /*
    ## Gradient of Es
    edge = node1 - node0

    edgeLen = np.linalg.norm(edge)
    tangent = edge / edgeLen
    epsX = edgeLen / l_k - 1
    dF_unit = EA * tangent * epsX
    dF = np.zeros((6))
    dF[0:3] = - dF_unit
    dF[3:6] = dF_unit

    ## Hessian of Es
    Id3 = np.eye(3)
    M = EA * ((1 / l_k - 1 / edgeLen) * Id3 + 1 / edgeLen * ( np.outer( edge, edge ) ) / edgeLen ** 2)

    dJ = np.zeros((6,6))
    dJ[0:3,0:3] = M
    dJ[3:6,3:6] = M
    dJ[0:3,3:6] = - M
    dJ[3:6,0:3] = - M
    return dF,dJ
    */

   // Gradient of Es
    Eigen::Vector3d edge = node1 - node0;
    double edgeLen = edge.norm();
    Eigen::Vector3d tangent = edge / edgeLen;
    double epsX = edgeLen / lk - 1.0;
    Eigen::Vector3d dF_unit = EA * tangent * epsX;
    
    Eigen::Map<
        Eigen::Matrix<
            double, 
            6, 1>
    > dFMap(dF);
    dFMap << -dF_unit, dF_unit;
   // 
   // // Hessian of Es
    
    Eigen::Matrix3d Id3 = Eigen::Matrix3d::Identity();

    // Note i did a transpose not in the original formula here...
    Eigen::Matrix3d M = EA * ((1 / lk - 1 / edgeLen) * Id3 + 1 / edgeLen * (edge * edge.transpose()) / (edgeLen * edgeLen));
    // Eigen::Vector3d dF = Eigen::Vector3d::Zero();
    Eigen::Map<Eigen::Matrix<double, 6, 6>> dJMap(dJ);
    dJMap.setZero();
    dJMap.block<3, 3>(0, 0) = M;
    dJMap.block<3, 3>(3, 3) = M;
    dJMap.block<3, 3>(0, 3) = -M;
    dJMap.block<3, 3>(3, 0) = -M;
    
}

static void calculateFs_Js(double* stretchingForcesVecPtr, SparseEntries* sparseEntries,// double* stretchingHessiansVecPtr, 
            const Eigen::Ref<const Eigen::Matrix<double, Eigen::Dynamic, 1>>& q, 
            const std::vector<glm::ivec2>& edgeVec, 
            const std::vector<double>& undeformedEdgeLengthVec, 
            const std::vector<double>& elasticStretchingVec, 
            uint32_t dofCount, uint32_t edgeCount,
            const std::span<uint32_t> oldToNewIndiceMapping) {

    Eigen::Map<
        Eigen::Matrix<
            double,
            Eigen::Dynamic, 1>
    > stretchForces(stretchingForcesVecPtr, dofCount, 1);

    // Eigen::Map<
    //     Eigen::Matrix<
    //     double,
    //     Eigen::Dynamic, Eigen::Dynamic
    //     >> stretchHessians(stretchingHessiansVecPtr, dofCount, dofCount);


    for (uint32_t edgeIdx = 0; edgeIdx < edgeCount; edgeIdx++) {
        glm::ivec2 kEdge = edgeVec[edgeIdx];
        uint32_t node0 = kEdge[0];
        uint32_t node1 = kEdge[1];
        double lk = undeformedEdgeLengthVec[edgeIdx];
        double ks  = elasticStretchingVec[edgeIdx];

        Eigen::Vector3d x0 = q.segment<3>(3 * node0);
        Eigen::Vector3d x1 = q.segment<3>(3 * node1);

        uint32_t ind[] = { 3 * node0, 3 * node0 + 1, 3 * node0 + 2,
                      3 * node1, 3 * node1 + 1, 3 * node1 + 2 };
        
        double dF[6];
        double dJ[36];
        Eigen::Map<
            Eigen::Matrix<
                double, 
                6, 1>
        > dFMap(dF);
        Eigen::Map<
            Eigen::Matrix<
                double, 
                6, 6>
        > dJMap(dJ);

        calculateGradEs_HessEs(x0, x1, lk, ks, dF, dJ);
        stretchForces(ind) = stretchForces(ind).eval() - dFMap;
        subFromSparseEntries(dJMap, ind, sparseEntries, oldToNewIndiceMapping);
        // stretchHessians(ind, ind) = stretchHessians(ind, ind).eval() - dJMap;

    }
}

static void SolveSystemWithConjugateGradient(const Eigen::Ref<const Eigen::SparseMatrix<double>>& J_free,
                        const Eigen::Ref<const Eigen::Matrix<double, Eigen::Dynamic, 1>>& f_free,
                        Eigen::Ref<Eigen::Matrix<double, Eigen::Dynamic, 1>> dq_free) {
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> solver;
    solver.compute(J_free);
    dq_free = solver.solve(f_free);
}

void ClothSystem::InitializePardiso(){
    pardisoData.mtype = 11;        /* Real symmetric? */

    // void* pt[64];

    // /* Pardiso control parameters. */
    // int      iparm[65];
    // double   dparm[64];
    // int      solver;
    // int      maxfct, mnum, phase, error, msglvl;

    // /* Number of processors. */
    // int      num_procs;

    // /* Auxiliary variables. */
    // char* var;
    // int      i;

    // double   ddum;              /* Double dummy */
    // int      idum;              /* Integer dummy. */

    /* -------------------------------------------------------------------- */
    /* ..  Setup Pardiso control parameters and initialize the solvers      */
    /*     internal adress pointers. This is only necessary for the FIRST   */
    /*     call of the PARDISO solver.                                      */
    /* ---------------------------------------------------------------------*/

    pardisoData.error = 0;
    pardisoData.solver = 0; /* use sparse direct solver */
    pardisoinit_d(pardisoData.pt, &pardisoData.mtype, &pardisoData.solver, &pardisoData.iparm[1], pardisoData.dparm, &pardisoData.error);

    pardisoData.num_procs = 12;
    pardisoData.iparm[3] = pardisoData.num_procs;
    pardisoData.iparm[7] = 1;
    pardisoData.iparm[11] = 0;
    pardisoData.iparm[13] = 0;
    pardisoData.iparm[27] = 0; /* recaling */

    pardisoData.maxfct = 1;         /* Maximum number of numerical factorizations.  */
    pardisoData.mnum = 1;         /* Which factorization to use. */

    pardisoData.msglvl = 0; // 0         /* Print statistical information  */
    pardisoData.error = 0;         /* Initialize error flag */
            /* Analyse and factorize the matrix */
    pardisoData.n = freeDOFCount;         /* Number of equations */
    pardisoData.nrhs = 1;         /* Number of right hand sides. */
}

static void SolveSystemWithPardiso(Eigen::Ref< Eigen::SparseMatrix<double>> J_free,
    Eigen::Ref<Eigen::Matrix<double, Eigen::Dynamic, 1>> f_free,
    Eigen::Ref<Eigen::Matrix<double, Eigen::Dynamic, 1>> dq_free,
    ClothSystem::PardisoData& pardisoData) {

    // int    n = J_free.rows();
    int    *ia = J_free.outerIndexPtr();
    int    *ja = J_free.innerIndexPtr();

    double  *a = J_free.valuePtr();

    double  *b = f_free.data();
    double  *x = dq_free.data();

    int      nnz = ia[pardisoData.n];




    /* -------------------------------------------------------------------- */
    /* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
    /*     notation.                                                        */
    /* -------------------------------------------------------------------- */
    shiftIndices(pardisoData.n, nnz, ia, ja, 1);


    /* -------------------------------------------------------------------- */
    /*  .. pardiso_chk_matrix(...)                                          */
    /*     Checks the consistency of the given matrix.                      */
    /*     Use this functionality only for debugging purposes               */
    /* -------------------------------------------------------------------- */
    //pardiso_chkmatrix_d(&pardisoData.mtype, &pardisoData.n, a, ia, ja, &pardisoData.error);
    //if (pardisoData.error != 0)
    //{
    //    std::printf("\nERROR in consistency of matrix: %d", pardisoData.error);
    //}


    // /* -------------------------------------------------------------------- */    
    // /* ..  Reordering and Symbolic Factorization.  This step also allocates */
    // /*     all memory that is necessary for the factorization.              */
    // /* -------------------------------------------------------------------- */ 
    
    double   ddum;              /* Double dummy */
    int      idum;              /* Integer dummy. */
    pardisoData.phase = 13;


    pardiso_d(pardisoData.pt, &pardisoData.maxfct, &pardisoData.mnum, &pardisoData.mtype, &pardisoData.phase,
        &pardisoData.n, a, ia, ja, &idum, &pardisoData.nrhs,
        &pardisoData.iparm[1], &pardisoData.msglvl, b, x, &pardisoData.error, pardisoData.dparm);

    if (pardisoData.error != 0)
    {
        std::printf("\nERROR during solve: %d", pardisoData.error);
        return;
    }



    /* -------------------------------------------------------------------- */
    /* ..  Convert matrix from 1-based Fortan notation to 0-based C         */
    /*     notation.                                                        */
    /* -------------------------------------------------------------------- */
    shiftIndices(pardisoData.n, nnz, ia, ja, -1);

    /* -------------------------------------------------------------------- */
    /* ..  Termination and release of memory.                               */
    /* -------------------------------------------------------------------- */
    pardisoData.phase = -1;                 /* Release internal memory. */

    pardiso_d(pardisoData.pt, &pardisoData.maxfct, &pardisoData.mnum, &pardisoData.mtype, &pardisoData.phase,
        &pardisoData.n, a, ia, ja, &idum, &pardisoData.nrhs,
        &pardisoData.iparm[1], &pardisoData.msglvl, &ddum, &ddum, &pardisoData.error, pardisoData.dparm);




    //std::printf("EXIT: Completed\n");s


}



void ClothSystem::CalculateForces() {

    auto currTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = currTime - lastTime;

    
    if (duration.count() > dt ) {
        Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> bufferStrides(4, 1);
        Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> tmpVecStrides(3, 1);
        // Iterate through every hinge
        Eigen::Map<
            Eigen::Matrix<
                double, 
                Eigen::Dynamic, 1>
        > massVec(massVector.data(), dofCount, 1);
        
        // Eigen::Map<Eigen::Matrix<
        //     double,
        //     Eigen::Dynamic, Eigen::Dynamic
        //     >> massMat(massMatrix.data(), dofCount, dofCount);
        // std::cout << massMat << std::endl;

        uint32_t iter = 0;
        double error = 10 * tol;

        Eigen::Map<
            Eigen::Matrix<
                double, 
                Eigen::Dynamic, 1>
        > externalForces(externalForcesVec.data(), dofCount, 1);


        // Eigen::MatrixXd positionX1 = positionMap.reshaped<Eigen::RowMajor>(dofCount, 1);

        // lastdofPositions = dofPositions;

        Eigen::Map<
            Eigen::Matrix<
                double, 
                Eigen::Dynamic, 1>
        > q0(lastdofPositions.data(), dofCount, 1);
        Eigen::Map<
            Eigen::Matrix<
                double, 
                Eigen::Dynamic, 1>
        > q(dofPositions.data(), dofCount, 1);
        Eigen::Map<
            Eigen::Matrix<
                double, 
                Eigen::Dynamic, 1>
        > u(dofVelocities.data(), dofCount, 1);

        q0 = q.eval();

        Eigen::Matrix<double, Eigen::Dynamic, 1> dq_free(freeDOFCount);

        // dofPositions = positionX1;
        while (error > tol) {
            if (iter > 5) {
                break;
            }
            
            auto startForceTime = std::chrono::high_resolution_clock::now();
            std::vector<double> bendingForcesVec(dofCount, 0.0);

            std::vector<double> stretchingForcesVec(dofCount, 0.0);

            Eigen::Map<
                Eigen::Matrix<
                    double, 
                    Eigen::Dynamic, 1>
            > bendingForces(bendingForcesVec.data(), dofCount, 1);

            Eigen::Map<
                Eigen::Matrix<
                    double, 
                    Eigen::Dynamic, 1>
            > stretchForces(stretchingForcesVec.data(), dofCount, 1);

            // Eigen::SparseMatrix<double> bendingSPMat(dofCount, dofCount);
            // Eigen::SparseMatrix<double> stretchSPMat(dofCount, dofCount);
            Eigen::SparseMatrix<double> jacobianForces(freeDOFCount, freeDOFCount);
            SparseEntries hessianEntries;
            
            //Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>> positionX1(positionMap, dofCount, 1);

            //Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>> velocityX1(velocityMap, dofCount, 1);

            //Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>> forceX1(forceMap, dofCount, 1);

            
            // Bending Forces -------------------------------------------------------------
            calculateFb_Jb(bendingForcesVec.data(), &hessianEntries, hingeVec, q, bendingStiffness, thetaBarVec, dofCount, hingeCount, oldToNewIndiceMapping);
            
            // Stretching Forces -------------------------------------------------------------
            calculateFs_Js(stretchingForcesVec.data(), &hessianEntries, q, edgeVec, undeformedEdgeLengthVec, elasticStretchingVec, dofCount, edgeCount, oldToNewIndiceMapping);
           
            auto endForceTime = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> forceTime = endForceTime - startForceTime;
            totalForceTime += forceTime.count();
            
            auto startSolverTime = std::chrono::high_resolution_clock::now();



            jacobianForces.setFromTriplets(hessianEntries.begin(), hessianEntries.end());
            
            
            Eigen::Matrix<double, Eigen::Dynamic, 1> totalForces = bendingForces + stretchForces + externalForces;
            Eigen::Matrix<double, Eigen::Dynamic, 1> f = (massVec / dt).cwiseProduct((q - q0) / dt - u) - totalForces;
            Eigen::Matrix<double, Eigen::Dynamic, 1> f_free = f(freeIndices);

            // Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> jacobianForces = bendingHessians + stretchHessians;
            // Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> J = massMatrix / (dt * dt) - jacobianForces;
            // Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> J_free = J(freeIndices, freeIndices);

            Eigen::SparseMatrix<double> J_free = massMatrix / (dt * dt) - jacobianForces;
            //std::cout << "Jacobian-----------------------------\n" << jacobianForces << std::endl;

            //std::cout  << "Jfree-----------------------------\n"  << J_free << std::endl;
            //Eigen::SparseMatrix<double> J_free = J(freeIndices, freeIndices);

            // double* f_free_ptr = f_free.data();
            // double* J_free_ptr = J_free.data();
            // std::cout << J_free << std::endl;
            //Eigen::SparseMatrix<double> sp_f_free(freeDOFCount, 1);
            //Eigen::SparseMatrix<double> sp_J_free(freeDOFCount, freeDOFCount);
            //
            //sp_f_free.reserve(Eigen::VectorXi::Constant(1, freeDOFCount));
            //sp_J_free.reserve(Eigen::VectorXi::Constant(1, freeDOFCount));


            // Eigen::Matrix<double, Eigen::Dynamic, 1> dq_free = J_free.ldlt().solve(f_free);
            J_free.makeCompressed();
            // SolveSystemWithConjugateGradient(J_free, f_free, dq_free);
            SolveSystemWithPardiso(J_free, f_free, dq_free, pardisoData);
            auto endSolverTime = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> solverTime = endSolverTime - startSolverTime;      
            totalSolveTime += solverTime.count();  
            
            q(freeIndices) = q(freeIndices).eval() - dq_free;

    

            error = f_free.cwiseAbs().sum();
            iter += 1;
            //std::cout << error << " ? " << tol << std::endl;
            
            
           //error = 0.0;
        }

        //fmt::print("Time: {:.6f}\n[", totalTime);
        //uint32_t perline = 0;
        //for (uint32_t freeIdx : freeIndices) {
        //   fmt::print("{: .8e} ", q(freeIdx));
        //   perline += 1;
        //   if (perline == 4) {
        //       perline = 0;
        //       fmt::print("\n ");
        //   }
        //}
        //fmt::print("]\n");

        u = (q - q0) / dt;

        for (uint32_t particleIdx = 0; particleIdx < particleCount; particleIdx++) {
            // dofPositions[dof] = q(dof, 0);
            // dofVelocities[dof] = u(dof, 0);
            clothPositionVec[particleIdx * 4] = dofPositions[particleIdx * 3];
            clothPositionVec[particleIdx * 4 + 1] = dofPositions[particleIdx * 3 + 1];
            clothPositionVec[particleIdx * 4 + 2] = dofPositions[particleIdx * 3 + 2];
        }

        glNamedBufferSubData(positionBuffer, 0, sizeof(float) * 4 * particleCount, clothPositionVec.data());
        lastTime = currTime;
        totalTime += dt;
        
    }
    
    //std::chrono::duration<double, std::milli> timeSinceStart = std::chrono::high_resolution_clock::now() - startTime;
    //if (timerRunning && timeSinceStart.count() > 5000) {
    //    fmt::print("t = 5: {} {} {}", dofPositions[particleCount * 3 - 3], dofPositions[particleCount * 3 - 2], dofPositions[particleCount * 3 - 1]);
    //    timerRunning = false;
    //}
    
}

void ClothSystem::BindBuffers() const{
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, PARTICLE_POSITIONS_SSBO_BINDING, positionBuffer);
    //glBindBufferBase(GL_SHADER_STORAGE_BUFFER, PARTICLE_VELOCITIES_SSBO_BINDING, velocityBuffer);
    //glBindBufferBase(GL_SHADER_STORAGE_BUFFER, PARTICLE_FORCES_SSBO_BINDING, forcesBuffer);
    //glBindBufferBase(GL_SHADER_STORAGE_BUFFER, PARTICLE_DATA_SSBO_BINDING, particleDataBuffer);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, PARTICLE_SYSTEM_DATA_SSBO_BINDING, particleSystemDataBuffer);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, PARTICLE_EDGE_DATA_SSBO_BINDING, edgesBuffer);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, PARTICLE_HINGE_DATA_SSBO_BINDING, hingesBuffer);
    glBindBufferBase(GL_UNIFORM_BUFFER, PICKED_DATA_UBO_BINDING, pickedDebugDataBuffer);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, PARTICLE_FIXED_NODES_SSBO_BINDING, fixedNodesDataBuffer);
    particleComputeShader->use();
}

void ClothSystem::SetupRender() {

    particleShader->use();
}


void ClothSystem::RenderDebug(GLuint VAO) {
    // glBindVertexArray(VAO);
    // glNamedBufferSubData(pickedDebugDataBuffer, 0, sizeof(PickedClothData), &pickedClothData);
    glNamedBufferSubData(fixedNodesDataBuffer, 0, sizeof(uint32_t)* fixedNodesCount, fixedNodes.data());

    //edgeDebugShader->use();
    //glVertexArrayVertexBuffer(VAO, 0, edgeDebugDataBuffer, 0, sizeof(glm::dvec4));
    //glVertexArrayElementBuffer(VAO, edgeDebugEBO);
    //glPointSize(10.0);
    //glLineWidth(10.0);
    //glDrawElements(GL_LINES, 2, GL_UNSIGNED_INT, 0);
    //glDrawElements(GL_POINTS, 2, GL_UNSIGNED_INT, 0);

    //hingeDebugShader->use();
    //glPointSize(20.0);
    //glLineWidth(20.0);
    //glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
    //glDrawElements(GL_POINTS, 6, GL_UNSIGNED_INT, 0);
    
    
    glBindVertexArray(fixedNodesVAO);
    fixedNodesDebugShader->use();
    glPointSize(20.0);
    glVertexArrayVertexBuffer(fixedNodesVAO, 0, fixedNodesDataBuffer, 0, sizeof(uint32_t));
    glDrawArrays(GL_POINTS, 0, fixedNodesCount);


    auto currTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = currTime - lastTime;
    
    //if (duration.count() > 120) {
    //    pickedClothData.pickedEdgeIdx = (pickedClothData.pickedEdgeIdx + 1) % edgeCount;
    //    // fmt::print("Edge {}: {} {}\n",
    //    //     pickedClothData.pickedEdgeIdx, 
    //    //     edgeVec[pickedClothData.pickedEdgeIdx].x, 
    //    //     edgeVec[pickedClothData.pickedEdgeIdx].y);

    //    // fmt::print("Hinge {}: {} {} {} {}\n",
    //    //     pickedClothData.pickedHingeIdx, 
    //    //     hingeVec[pickedClothData.pickedHingeIdx].x, 
    //    //     hingeVec[pickedClothData.pickedHingeIdx].y,
    //    //     hingeVec[pickedClothData.pickedHingeIdx].z,
    //    //     hingeVec[pickedClothData.pickedHingeIdx].w);
    //    pickedClothData.pickedHingeIdx = (pickedClothData.pickedHingeIdx + 1) % hingeCount;        
    //    lastTime = currTime;
    //}


}