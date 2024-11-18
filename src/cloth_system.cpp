#include <cloth_system.hpp>
#include <shader_c.hpp>
#include <shader_s.hpp>
#include <chrono>
#include <iostream>
#include <Eigen/Dense>
 

//ClothSystem::ClothSystem(ClothSystemParameters* params) :
//    BaseParticleSystem(params) {
//
//};

void ClothSystem::InitializeSystemData(void* params) {
    ClothSystemParameters* clothParams = static_cast<ClothSystemParameters*>(params);
    clothParams->particleCount = (clothParams->clothSideLength + 1) * (clothParams->clothSideLength + 1);
    particleSystemDataBlock.particleCount = clothParams->particleCount;
    particleCount = clothParams->particleCount;
    particleSystemDataBlock.timeStep = clothParams->timeStep;
    dofCount = 3 * particleCount;
    vecValCount = particleCount * 4;
    
    indiceCount = 6 * clothParams->clothSideLength * clothParams->clothSideLength;
    //particleSystemDataBlock.tolerance = clothParams->tolerance;
    //particleSystemDataBlock.fluidViscosity = clothParams->fluidViscosity;
    //particleSystemDataBlock.gravityAccel = clothParams->gravityAccel;
    //particleSystemDataBlock.youngModulus = 1.0e7f;
    //particleSystemDataBlock.thickness = 1.0e3f;
    //particleSystemDataBlock.bendingStiffness = 2.0f / sqrt(3.0f) * particleSystemDataBlock.youngModulus * pow(particleSystemDataBlock.thickness, 3.0f) / 12.0f;
    //particleSystemDataBlock.stretchStiffness = 1.0f;
};



void ClothSystem::InitializeBufferData(void* params) {
    ClothSystemParameters* clothParams = static_cast<ClothSystemParameters*>(params);
    clothParams->particleCount = (clothParams->clothSideLength+ 1) * (clothParams->clothSideLength + 1);
    // count particles
    positionVec.resize(vecValCount);
    velocityVec.resize(vecValCount);
    forceVec.resize(vecValCount);
    normalsVec.resize(clothParams->particleCount);
    particleDataVec.resize(clothParams->particleCount);

    // tangentsVec.resize(clothParams->particleCount);
    // bitangentsVec.resize(clothParams->particleCount);


    glm::vec4 meshStartPosition = glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);
    float yIncrement = 0.5f * sqrt(3.0f) * clothParams->cellSideLength;


    uint32_t clothSideParticleCount = clothParams->clothSideLength + 1;
    float deltaMass = clothParams->totalMass / clothParams->particleCount;

    for (uint32_t rowIdx = 0; rowIdx < clothSideParticleCount; ++rowIdx) {
        glm::vec4 rowStartPosition =
            meshStartPosition;

        rowStartPosition.z -= yIncrement * rowIdx;
        rowStartPosition.y = 0.5f*sin(rowIdx * clothParams->cellSideLength * 2.0f);
        if (rowIdx % 2 == 1) {
            rowStartPosition.x += 0.5f * clothParams->cellSideLength;
        }

        glm::vec4 currentPosition = rowStartPosition;
        for (uint32_t colIdx = 0; colIdx < clothSideParticleCount; ++colIdx) {
            currentPosition.y = rowStartPosition.y + 0.5f * sin(colIdx * clothParams->cellSideLength * 2.0f);
            uint32_t particleIdx = rowIdx * clothSideParticleCount + colIdx;
            for (uint32_t i = 0; i < 4; ++i) {
                positionVec[particleIdx * 4 + i] = currentPosition[i];
                velocityVec[particleIdx * 4 + i] = 0.0f;
                forceVec[particleIdx * 4 + i] = 0.0f;
            }
            // positionVec[particleIdx] = currentPosition;
            // velocityVec[particleIdx] = glm::vec4(0.0f);
            // forceVec[particleIdx] = glm::vec4(0.0f);
            particleDataVec[particleIdx].mass = deltaMass;

            // normalsVec[particleIdx] = glm::vec4(0.0f);
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

                edgeVec[edgeIdx + 0] = glm::vec2(aIdx, dIdx);
                edgeVec[edgeIdx + 1] = glm::vec2(dIdx, cIdx);
                edgeVec[edgeIdx + 2] = glm::vec2(cIdx, aIdx);
                
            }
            else {
                indicesVec[indiceIdx + 0] = aIdx;
                indicesVec[indiceIdx + 1] = bIdx;
                indicesVec[indiceIdx + 2] = cIdx;

                indicesVec[indiceIdx + 3] = bIdx;
                indicesVec[indiceIdx + 4] = dIdx;
                indicesVec[indiceIdx + 5] = cIdx;
                
                hingeVec[hingeIdx] = glm::ivec4(cIdx, bIdx, dIdx, aIdx);
            
                edgeVec[edgeIdx + 0] = glm::vec2(bIdx, dIdx);
                edgeVec[edgeIdx + 1] = glm::vec2(dIdx, cIdx);
                edgeVec[edgeIdx + 2] = glm::vec2(cIdx, bIdx);

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




            /*fmt::print("\nFirst Triangle\n");
            for (uint32_t j = 0; j < 3; ++j) {
                fmt::print("{} ", indicesVec[indiceIdx + j]);
            }
            fmt::print("\nSecond Triangle\n");
            for (uint32_t j = 3; j < 6; ++j) {
                fmt::print("{} ", indicesVec[indiceIdx + j]);
            }*/
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
    const float Y = 1.0e7f;
    const float h = 1.0e3f;
    bendingStiffness = 2.0f / sqrt(3.0f) * Y * pow(h, 3.0f) / 12.0f;
    for (uint32_t edgeIdx = 0; edgeIdx < edgeCount; ++edgeIdx) {
        glm::ivec2 edge = edgeVec[edgeIdx];
        glm::vec4 edgeStart (positionVec[4 * edge[0]], positionVec[4 * edge[0] + 1], positionVec[4 * edge[0] + 2], 1.0f);
        glm::vec4 edgeEnd (positionVec[4 * edge[1]], positionVec[4 * edge[1] + 1], positionVec[4 * edge[1] + 2], 1.0f);
        glm::vec4 edgeVec = edgeEnd - edgeStart;
        undeformedEdgeLengthVec[edgeIdx] = glm::length(edgeVec);
        elasticStretchingVec[edgeIdx] = 0.5f * sqrt(3.0) * Y * h * pow(undeformedEdgeLengthVec[edgeIdx], 2.0f);
        
    }
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

    particleShader = new Shader(vertPath.c_str(), fragPath.c_str());
    particleComputeShader = new ComputeShader(compPath.c_str());

    edgeDebugShader = new Shader(edgeDebugVertPath.c_str(), edgeDebugFragPath.c_str());
    hingeDebugShader = new Shader(hingeDebugVertPath.c_str(), hingeDebugFragPath.c_str());

}

void ClothSystem::InitializeBuffers(){
    BaseParticleSystem::InitializeBuffers();
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
       sizeof(glm::vec4) * 4, 
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

    //https://stackoverflow.com/questions/12399422/how-to-set-linker-flags-for-openmp-in-cmakes-try-compile-function
    lastTime = std::chrono::high_resolution_clock::now();
  Eigen::MatrixXd m(2,2);
  m(0,0) = 3;
  m(1,0) = 2.5;
  m(0,1) = -1;
  m(1,1) = m(1,0) + m(0,1);
  //fmt::print("Eigen used here: {}", m);
  std::cout << m << std::endl;


}

    



/*
"""# Bending energy for a shell, it's gradient, and Hessian"""


def getEb_Shell(x0, x1=None, x2=None, x3=None, theta_bar=0, kb=1.0):
    """
    Compute the bending energy for a shell.

    Returns:
    E (scalar): Bending energy.
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

    # E = 0.5 * kb * (theta-thetaBar)^2
    E = 0.5 * kb * (theta - theta_bar) ** 2

    return E
    */

    /*

def objfun(qGuess, q0, u, freeIndex, dt, tol,
                massVector, mMat,
                ks, lk, edges,
                kb, thetaBar, hinges,
                Fg):

  q = qGuess
  ndof = len(q)
  iter = 0
  error = 10 * tol

  while error > tol:

    # Calculating Bending
    Fb = np.zeros(ndof)
    Jb = np.zeros((ndof,ndof))
    for kHinge in range(hinges.shape[0]):
      node0 = hinges[kHinge,0]
      node1 = hinges[kHinge,1]
      node2 = hinges[kHinge,2]
      node3 = hinges[kHinge,3]
      x0 = q[3*node0:3*node0+3]
      x1 = q[3*node1:3*node1+3]
      x2 = q[3*node2:3*node2+3]
      x3 = q[3*node3:3*node3+3]
      ind = [3*node0, 3*node0 + 1, 3*node0 + 2,
             3*node1, 3*node1 + 1, 3*node1 + 2,
             3*node2, 3*node2 + 1, 3*node2 + 2,
             3*node3, 3*node3 + 1, 3*node3 + 2]
      dF, dJ = gradEb_hessEb_Shell(x0, x1, x2, x3, thetaBar, kb)
      Fb[ind] = Fb[ind] - dF
      Jb[np.ix_(ind,ind)] = Jb[np.ix_(ind,ind)] - dJ

    # Calculating stretching
    Fs = np.zeros(ndof)
    Js = np.zeros((ndof,ndof))
    for kEdge in range(edges.shape[0]):
      node0 = edges[kEdge,0]
      node1 = edges[kEdge,1]
      x0 = q[3*node0:3*node0+3]
      x1 = q[3*node1:3*node1+3]
      ind = [3*node0, 3*node0 + 1, 3*node0 + 2,
             3*node1, 3*node1 + 1, 3*node1 + 2]
      dF, dJ = gradEs_hessEs(x0, x1, lk[kEdge], ks[kEdge])
      Fs[ind] = Fs[ind] - dF
      Js[np.ix_(ind,ind)] = Js[np.ix_(ind,ind)] - dJ

    # Calculating total force
    F = Fg + Fb + Fs # Viscous forces can sometimes be useful
    JForces = Js + Jb

    f = massVector/dt * ( (q-q0)/dt - u ) - F
    J = mMat / dt ** 2 - JForces

    f_free = f[freeIndex]
    J_free = J[np.ix_(freeIndex, freeIndex)]
    dq_free = np.linalg.solve(J_free, f_free)

    q[freeIndex] = q[freeIndex] - dq_free

    error = np.sum( np.abs(f_free) )
    iter += 1

    print('Iter = %d' % iter)
    print('Error = %f' % error)

  u = (q - q0) / dt

  return q, u


*/


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
    #   angle: float, the signed angle (in radians) from vector "u" to vector "v".
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

static float signedAngle(const Eigen::RowVector3f& u, const Eigen::RowVector3f& v, const Eigen::RowVector3f& n){
    Eigen::RowVector3f w = u.cross(v);
    float angle = atan2(w.norm(), u.dot(v));
    if (n.dot(w) < 0){
        angle = -angle;
    }
    return angle;
}


static Eigen::Matrix3f mmt(const Eigen::Matrix3f& matrix){
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

static float getTheta(const Eigen::RowVector3f& node0, const Eigen::RowVector3f& node1, const Eigen::RowVector3f& node2, const Eigen::RowVector3f& node3){
    Eigen::RowVector3f m_e0 = node1 - node0;
    Eigen::RowVector3f m_e1 = node2 - node0;
    Eigen::RowVector3f m_e2 = node3 - node0;

    Eigen::RowVector3f n0 = m_e0.cross(m_e1);
    Eigen::RowVector3f n1 = m_e2.cross(m_e0);

    return signedAngle(n0, n1, m_e0);
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

static Eigen::MatrixXf uvT(const Eigen::RowVector3f& u, const Eigen::RowVector3f& v){
    return u.transpose() * v;
}

static void hessTheta(const Eigen::RowVector3f& node0, const Eigen::RowVector3f& node1, const Eigen::RowVector3f& node2, const Eigen::RowVector3f& node3, float* hess){
    Eigen::RowVector3f m_e0 = node1 - node0;
    Eigen::RowVector3f m_e1 = node2 - node0;
    Eigen::RowVector3f m_e2 = node3 - node0;
    Eigen::RowVector3f m_e3 = node2 - node1;
    Eigen::RowVector3f m_e4 = node3 - node1;

    float m_cosA1 = m_e0.dot(m_e1) / (m_e0.norm() * m_e1.norm());
    float m_cosA2 = m_e0.dot(m_e2) / (m_e0.norm() * m_e2.norm());
    float m_cosA3 = -m_e0.dot(m_e3) / (m_e0.norm() * m_e3.norm());
    float m_cosA4 = -m_e0.dot(m_e4) / (m_e0.norm() * m_e4.norm());

    float m_sinA1 = m_e0.cross(m_e1).norm() / (m_e0.norm() * m_e1.norm());
    float m_sinA2 = m_e0.cross(m_e2).norm() / (m_e0.norm() * m_e2.norm());
    float m_sinA3 = -m_e0.cross(m_e3).norm() / (m_e0.norm() * m_e3.norm());
    float m_sinA4 = -m_e0.cross(m_e4).norm() / (m_e0.norm() * m_e4.norm());

    Eigen::RowVector3f m_nn1 = m_e0.cross(m_e3);
    m_nn1 = m_nn1 / m_nn1.norm();
    Eigen::RowVector3f m_nn2 = -m_e0.cross(m_e4);
    m_nn2 = m_nn2 / m_nn2.norm();

    float m_h1 = m_e0.norm() * m_sinA1;
    float m_h2 = m_e0.norm() * m_sinA2;
    float m_h3 = -m_e0.norm() * m_sinA3;
    float m_h4 = -m_e0.norm() * m_sinA4;
    float m_h01 = m_e1.norm() * m_sinA1;
    float m_h02 = m_e2.norm() * m_sinA2;

    Eigen::RowVector3f m_m1 = m_nn1.cross(m_e1) / m_e1.norm();
    Eigen::RowVector3f m_m2 = -m_nn2.cross(m_e2) / m_e2.norm();
    Eigen::RowVector3f m_m3 = -m_nn1.cross(m_e3) / m_e3.norm();
    Eigen::RowVector3f m_m4 = m_nn2.cross(m_e4) / m_e4.norm();
    Eigen::RowVector3f m_m01 = -m_nn1.cross(m_e0) / m_e0.norm();
    Eigen::RowVector3f m_m02 = m_nn2.cross(m_e0) / m_e0.norm();
    //
    Eigen::Map<Eigen::MatrixXf> hess_theta(hess, 12, 12);
    hess_theta.setZero();

    Eigen::MatrixXf M331 = m_cosA3 / (m_h3 * m_h3) * uvT(m_m3, m_nn1);
    Eigen::MatrixXf M311 = m_cosA3 / (m_h3 * m_h1) * uvT(m_m1, m_nn1);
    Eigen::MatrixXf M131 = m_cosA1 / (m_h1 * m_h3) * uvT(m_m3, m_nn1);
    Eigen::MatrixXf M3011 = m_cosA3 / (m_h3 * m_h01) * uvT(m_m01, m_nn1);
    Eigen::MatrixXf M111 = m_cosA1 / (m_h1 * m_h1) * uvT(m_m1, m_nn1);
    Eigen::MatrixXf M1011 = m_cosA1 / (m_h1 * m_h01) * uvT(m_m01, m_nn1);
    
    Eigen::MatrixXf M442 = m_cosA4 / (m_h4 * m_h4) * uvT(m_m4, m_nn2);
    Eigen::MatrixXf M422 = m_cosA4 / (m_h4 * m_h2) * uvT(m_m2, m_nn2);
    Eigen::MatrixXf M242 = m_cosA2 / (m_h2 * m_h4) * uvT(m_m4, m_nn2);
    Eigen::MatrixXf M4022 = m_cosA4 / (m_h4 * m_h02) * uvT(m_m02, m_nn2);
    Eigen::MatrixXf M222 = m_cosA2 / (m_h2 * m_h2) * uvT(m_m2, m_nn2);
    Eigen::MatrixXf M2022 = m_cosA2 / (m_h2 * m_h02) * uvT(m_m02, m_nn2);

    Eigen::MatrixXf B1 = 1 / m_e0.norm() / m_e0.norm() * uvT(m_nn1, m_m01);
    Eigen::MatrixXf B2 = 1 / m_e0.norm() / m_e0.norm() * uvT(m_nn2, m_m02);

    Eigen::MatrixXf N13 = 1 / (m_h01 * m_h3) * uvT(m_nn1, m_m3);
    Eigen::MatrixXf N24 = 1 / (m_h02 * m_h4) * uvT(m_nn2, m_m4);
    Eigen::MatrixXf N11 = 1 / (m_h01 * m_h1) * uvT(m_nn1, m_m1);
    Eigen::MatrixXf N22 = 1 / (m_h02 * m_h2) * uvT(m_nn2, m_m2);
    Eigen::MatrixXf N101 = 1 / (m_h01 * m_h01) * uvT(m_nn1, m_m01);
    Eigen::MatrixXf N202 = 1 / (m_h02 * m_h02) * uvT(m_nn2, m_m02);

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

static void gradTheta(const Eigen::RowVector3f& node0, const Eigen::RowVector3f& node1, const Eigen::RowVector3f& node2, const Eigen::RowVector3f& node3, float* grad){
    Eigen::RowVector3f m_e0 = node1 - node0;
    Eigen::RowVector3f m_e1 = node2 - node0;
    Eigen::RowVector3f m_e2 = node3 - node0;
    Eigen::RowVector3f m_e3 = node2 - node1;
    Eigen::RowVector3f m_e4 = node3 - node1;

    float m_cosA1 = m_e0.dot(m_e1) / (m_e0.norm() * m_e1.norm());
    float m_cosA2 = m_e0.dot(m_e2) / (m_e0.norm() * m_e2.norm());
    float m_cosA3 = -m_e0.dot(m_e3) / (m_e0.norm() * m_e3.norm());
    float m_cosA4 = -m_e0.dot(m_e4) / (m_e0.norm() * m_e4.norm());

    float m_sinA1 = m_e0.cross(m_e1).norm() / (m_e0.norm() * m_e1.norm());
    float m_sinA2 = m_e0.cross(m_e2).norm() / (m_e0.norm() * m_e2.norm());
    float m_sinA3 = -m_e0.cross(m_e3).norm() / (m_e0.norm() * m_e3.norm());
    float m_sinA4 = -m_e0.cross(m_e4).norm() / (m_e0.norm() * m_e4.norm());

    Eigen::RowVector3f m_nn1 = m_e0.cross(m_e3);
    m_nn1 /= m_nn1.norm();
    Eigen::RowVector3f m_nn2 = -m_e0.cross(m_e4);
    m_nn2 /= m_nn2.norm();

    float m_h1 = m_e0.norm() * m_sinA1;
    float m_h2 = m_e0.norm() * m_sinA2;
    float m_h3 = -m_e0.norm() * m_sinA3;
    float m_h4 = -m_e0.norm() * m_sinA4;
    float m_h01 = m_e1.norm() * m_sinA1;
    float m_h02 = m_e2.norm() * m_sinA2;

    Eigen::RowVector3f gradTheta0 = m_cosA3 * m_nn1 / m_h3 + m_cosA4 * m_nn2 / m_h4;
    Eigen::RowVector3f gradTheta1 = m_cosA1 * m_nn1 / m_h1 + m_cosA2 * m_nn2 / m_h2;
    Eigen::RowVector3f gradTheta2 = -m_nn1 / m_h01;
    Eigen::RowVector3f gradTheta3 = -m_nn2 / m_h02;

    Eigen::Map<Eigen::RowVectorXf> dFMap(grad, 12);
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


static void calculateGradEb_HessEb_Shell(const Eigen::RowVector3f& node0, const Eigen::RowVector3f& node1, const Eigen::RowVector3f& node2, const Eigen::RowVector3f& node3, float thetaBar, float kb, float* dF, float* dJ){
    /*
def gradEb_hessEb_Shell(x0, x1=None, x2=None, x3=None, theta_bar=0, kb=1.0):
    """
    Compute the gradient and Hessian of the bending energy for a shell.

    Parameters:
    x0 (array): Can either be a 3-element array (single point) or a 12-element array.
    x1, x2, x3 (arrays): Optional, 3-element arrays specifying points.
    theta_bar (float): Reference angle.
    kb (float): Bending stiffness.

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
    float theta = getTheta(node0, node1, node2, node3);
    float grad[12] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
    gradTheta(node0, node1, node2, node3, grad);

    Eigen::Map<Eigen::RowVectorXf> gradMap(grad, 12);
    Eigen::Map<Eigen::RowVectorXf> dFMap(dF, 12);

    dFMap = 0.5 * kb * (2 * (theta - thetaBar) * gradMap);
    
    float hess[144];
    hessTheta(node0, node1, node2, node3, hess);
    Eigen::Map<Eigen::Matrix<float, 12, 12, Eigen::RowMajor>> dJMap(dJ);
    Eigen::Map<Eigen::Matrix<float, 12, 12, Eigen::RowMajor>> hessMap(hess);
    //// TODO Which is the correct transpose?
    dJMap = 0.5 * kb * (2 *  gradMap.transpose() * gradMap + 2 * (theta - thetaBar) * hessMap);
}

static void calculateFb_Jb(float* bendingForcesVecPtr, float* bendingHessiansVecPtr, const std::vector<glm::ivec4>& hingeVec, const Eigen::MatrixXf& positionMap, float kb, float thetaBar, uint32_t dofCount, uint32_t hingeCount) {
    Eigen::Map<
        Eigen::Matrix<
        float,
        Eigen::Dynamic, Eigen::Dynamic,
        Eigen::RowMajor
        >, Eigen::Unaligned> bendingForces(bendingForcesVecPtr, dofCount, 1);
    Eigen::Map<
        Eigen::Matrix<
        float,
        Eigen::Dynamic, Eigen::Dynamic,
        Eigen::RowMajor
        >, Eigen::Unaligned> bendingHessians(bendingHessiansVecPtr, dofCount, dofCount);


    // Bending Forces
    for (uint32_t hingeIdx = 0; hingeIdx < hingeCount; ++hingeIdx) {
        glm::ivec4 kHinge = hingeVec[hingeIdx];
        uint32_t node0 = kHinge[0];
        uint32_t node1 = kHinge[1];
        uint32_t node2 = kHinge[2];
        uint32_t node3 = kHinge[3];

        Eigen::RowVector3f x0 = positionMap.row(node0);
        Eigen::RowVector3f x1 = positionMap.row(node1);
        Eigen::RowVector3f x2 = positionMap.row(node2);
        Eigen::RowVector3f x3 = positionMap.row(node3);

        int ind[] = { 3 * node0, 3 * node0 + 1, 3 * node0 + 2,
                      3 * node1, 3 * node1 + 1, 3 * node1 + 2,
                      3 * node2, 3 * node2 + 1, 3 * node2 + 2,
                      3 * node3, 3 * node3 + 1, 3 * node3 + 2 };

        float df[12] = { -1.0f, -2.0f, -3.0f, 1.0f, 2.0f, 3.0f, -4.0f, -5.0f, -6.0f, 4.0f, 5.0f, 6.0f };
        float dj[144];
        for (int i = 0; i < 144; i++) {
            dj[i] = i;
        }

        calculateGradEb_HessEb_Shell(x0, x1, x2, x3, 0.0f, 1.0f, df, dj);

        bendingForces(ind, Eigen::placeholders::all) = bendingForces(ind, Eigen::placeholders::all) - Eigen::Map<Eigen::Matrix<float, 12, 1>>(df);
        bendingHessians(ind, ind) = bendingHessians(ind, ind) - Eigen::Map<Eigen::Matrix<float, 12, 12>>(dj);

        //    //fmt
    }
}


static void calculateGradEs_HessEs(const Eigen::RowVector3f& node0, const Eigen::RowVector3f& node1, float lk, float EA, float* dF, float* dJ) {
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
    Eigen::RowVector3f edge = node1 - node0;
    float edgeLen = edge.norm();
    Eigen::RowVector3f tangent = edge / edgeLen;
    float epsX = edgeLen / lk - 1.0f;
    Eigen::RowVector3f dF_unit = EA * tangent * epsX;
    Eigen::Map<Eigen::RowVectorXf> dFMap(dF, 6);
    dFMap << -dF_unit, dF_unit;
   // 
   // // Hessian of Es
    Eigen::Matrix3f Id3 = Eigen::Matrix3f::Identity();

    // Note i did a transpose not in the original formula here...
    Eigen::Matrix3f M = EA * ((1 / lk - 1 / edgeLen) * Id3 + 1 / edgeLen * (edge.transpose() * edge) / (edgeLen * edgeLen));
    // Eigen::RowVector3f dF = Eigen::RowVector3f::Zero();
    Eigen::Map<Eigen::Matrix<float, 6, 6, Eigen::RowMajor>> dJMap(dJ);
    dJMap.setZero();
    dJMap.block<3, 3>(0, 0) = M;
    dJMap.block<3, 3>(3, 3) = M;
    dJMap.block<3, 3>(0, 3) = -M;
    dJMap.block<3, 3>(3, 0) = -M;
}

static void calculateFs_Js(float* stretchingForcesVecPtr, float* stretchingHessiansVecPtr, const Eigen::MatrixXf& positions, const std::vector<glm::ivec2>& edgeVec, const std::vector<float>& undeformedEdgeLengthVec, const std::vector<float>& elasticStretchingVec, uint32_t dofCount, uint32_t edgeCount) {

    Eigen::Map<
        Eigen::Matrix<
        float,
        Eigen::Dynamic, Eigen::Dynamic,
        Eigen::RowMajor
        >, Eigen::Unaligned> stretchForces(stretchingForcesVecPtr, dofCount, 1);

    Eigen::Map<
        Eigen::Matrix<
        float,
        Eigen::Dynamic, Eigen::Dynamic,
        Eigen::RowMajor
        >, Eigen::Unaligned> stretchHessians(stretchingHessiansVecPtr, dofCount, dofCount);


    for (uint32_t edgeIdx = 0; edgeIdx < edgeCount; edgeIdx++) {
        glm::ivec2 kEdge = edgeVec[edgeIdx];
        uint32_t node0 = kEdge[0];
        uint32_t node1 = kEdge[1];
        float lk = undeformedEdgeLengthVec[edgeIdx];
        float ks  = elasticStretchingVec[edgeIdx];

        Eigen::RowVector3f x0 = positions.row(node0);
        Eigen::RowVector3f x1 = positions.row(node1);

        int ind[] = { 3 * node0, 3 * node0 + 1, 3 * node0 + 2,
                      3 * node1, 3 * node1 + 1, 3 * node1 + 2 };

        float df[6] = {-1.0f, -2.0f, -3.0f, 1.0f, 2.0f, 3.0f};
        float dj[36];
        for (int i = 0; i < 36; i++) {
            dj[i] = i;
        }

        calculateGradEs_HessEs(x0, x1, lk, ks, df, dj);
        stretchForces(ind, Eigen::placeholders::all) = stretchForces(ind, Eigen::placeholders::all) -  Eigen::Map<Eigen::Matrix<float, 6, 1>>(df);
        stretchHessians(ind, ind) = stretchHessians(ind, ind) - Eigen::Map<Eigen::Matrix<float, 6, 6>>(dj);

        //fmt::print("Edge {} {}\n", node0, node1);
        //std::cout << "Vec0s ----------- " << x0 << std::endl;
        //std::cout << "Vec1s ----------- " << x1 << std::endl;

        //std::cout << "Stretching Forces -----------\n" << stretchForces << std::endl;
        //std::cout << "Stretching Hessians -----------\n" << stretchHessians << std::endl;
    }

}

void ClothSystem::CalculateForces() {
    // Eigen::Map<
    //     Eigen::MatrixXf, 
    //     Eigen::Unaligned, 
    //     Eigen::Stride<0, 0>    
    // > positionMap (positionVec.data(), 1, vecValCount, Eigen::Stride<0, 0>(0, 0));
    // //std::cout << "Vec--------\n" << positionMap << std::endl;

    // Eigen::Map<
    //     Eigen::MatrixXf,
    //     Eigen::Unaligned,
    //     Eigen::Stride<1, 4>
    // > testMap(positionVec.data(), 3, particleCount, Eigen::Stride<1, 4>(1, 4));
    // //std::cout << "Mat--------\n" << testMap << std::endl;
     
     
    // Define the float array (input data)
    // float rawArray[] = {
    //     0, 0, 0, 1,
    //     0.1, 0.0993347, 0, 1,
    //     0.05, 0.0993347, -0.0866025, 1,
    //     0.15, 0.198669, -0.0866025, 1
    // };



    // Define strides: 4 elements per row, skip 1 element per column
    Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> bufferStrides(4, 1);
    Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> tmpVecStrides(3, 1);


    Eigen::Map<
        Eigen::Matrix<
            float, 
            Eigen::Dynamic, 3, 
            Eigen::RowMajor
        >, Eigen::Unaligned, 
        Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>
    > positionMap(positionVec.data(), particleCount, 3, bufferStrides);


    // Iterate through every hinge
    //std::vector<float> bendingForces(dofCount, 0.0f);
    std::vector<float> bendingForcesVec(dofCount, 0.0f);
    std::vector<float> bendingHessiansVec(dofCount * dofCount, 0.0f);

    calculateFb_Jb(bendingForcesVec.data(), bendingHessiansVec.data(), hingeVec, positionMap, bendingStiffness, 0.0f, dofCount, hingeCount);

        // Fb = np.zeros(ndof)
        // Jb = np.zeros((ndof,ndof))
        // for kHinge in range(hinges.shape[0]):
        //   node0 = hinges[kHinge,0]
        //   node1 = hinges[kHinge,1]
        //   node2 = hinges[kHinge,2]
        //   node3 = hinges[kHinge,3]
        //   x0 = q[3*node0:3*node0+3]
        //   x1 = q[3*node1:3*node1+3]
        //   x2 = q[3*node2:3*node2+3]
        //   x3 = q[3*node3:3*node3+3]
        //   ind = [3*node0, 3*node0 + 1, 3*node0 + 2,
        //          3*node1, 3*node1 + 1, 3*node1 + 2,
        //          3*node2, 3*node2 + 1, 3*node2 + 2,
        //          3*node3, 3*node3 + 1, 3*node3 + 2]
        //   dF, dJ = gradEb_hessEb_Shell(x0, x1, x2, x3, thetaBar, kb)
        //   Fb[ind] = Fb[ind] - dF
        //   Jb[np.ix_(ind,ind)] = Jb[np.ix_(ind,ind)] - dJ
// Stretching Forces -------------------------------------------------------------

    //// Map the raw array to a 4x3 matrix with the defined bufferStrides
    //Eigen::Map<Eigen::Matrix<float, 4, 3, Eigen::RowMajor>, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>
    //    positions(rawArray, 4, 3, bufferStrides);

    // Print the extracted positions
    //std::cout << "Extracted 4x3 Positions Matrix:\n" << positions << "\n";
    // Iterate through every edge

    std::vector<float> stretchingForcesVec(dofCount, 0.0f);
    std::vector<float> stretchingHessiansVec(dofCount * dofCount, 0.0f);
    calculateFs_Js(stretchingForcesVec.data(), stretchingHessiansVec.data(), positionMap, edgeVec, undeformedEdgeLengthVec, elasticStretchingVec, dofCount, edgeCount);
    /*

  q = qGuess
  ndof = len(q)
  iter = 0
  error = 10 * tol

  while error > tol:

    # Calculating Bending
    Fb = np.zeros(ndof)
    Jb = np.zeros((ndof,ndof))
    for kHinge in range(hinges.shape[0]):
      node0 = hinges[kHinge,0]
      node1 = hinges[kHinge,1]
      node2 = hinges[kHinge,2]
      node3 = hinges[kHinge,3]
      x0 = q[3*node0:3*node0+3]
      x1 = q[3*node1:3*node1+3]
      x2 = q[3*node2:3*node2+3]
      x3 = q[3*node3:3*node3+3]
      ind = [3*node0, 3*node0 + 1, 3*node0 + 2,
             3*node1, 3*node1 + 1, 3*node1 + 2,
             3*node2, 3*node2 + 1, 3*node2 + 2,
             3*node3, 3*node3 + 1, 3*node3 + 2]
      dF, dJ = gradEb_hessEb_Shell(x0, x1, x2, x3, thetaBar, kb)
      Fb[ind] = Fb[ind] - dF
      Jb[np.ix_(ind,ind)] = Jb[np.ix_(ind,ind)] - dJ

    # Calculating stretching
    Fs = np.zeros(ndof)
    Js = np.zeros((ndof,ndof))
    for kEdge in range(edges.shape[0]):
      node0 = edges[kEdge,0]
      node1 = edges[kEdge,1]
      x0 = q[3*node0:3*node0+3]
      x1 = q[3*node1:3*node1+3]
      ind = [3*node0, 3*node0 + 1, 3*node0 + 2,
             3*node1, 3*node1 + 1, 3*node1 + 2]
      dF, dJ = gradEs_hessEs(x0, x1, lk[kEdge], ks[kEdge])
      Fs[ind] = Fs[ind] - dF
      Js[np.ix_(ind,ind)] = Js[np.ix_(ind,ind)] - dJ

    # Calculating total force
    F = Fg + Fb + Fs # Viscous forces can sometimes be useful
    JForces = Js + Jb

    f = massVector/dt * ( (q-q0)/dt - u ) - F
    J = mMat / dt ** 2 - JForces

    f_free = f[freeIndex]
    J_free = J[np.ix_(freeIndex, freeIndex)]
    dq_free = np.linalg.solve(J_free, f_free)

    q[freeIndex] = q[freeIndex] - dq_free

    error = np.sum( np.abs(f_free) )
    iter += 1

    print('Iter = %d' % iter)
    print('Error = %f' % error)

  u = (q - q0) / dt

  return q, u
    
    */


    // glNamedBufferSubData(
    //     particleSystemDataBuffer,
    //     0,
    //     sizeof(ParticleSystemDataBlock),
    //     &particleSystemDataBlock
    // );
    // glDispatchCompute(particleCount, 1u, 1u);
    // glMemoryBarrier(GL_ALL_BARRIER_BITS);

    // glFinish();


}

void ClothSystem::BindBuffers() const{
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, PARTICLE_POSITIONS_SSBO_BINDING, positionBuffer);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, PARTICLE_VELOCITIES_SSBO_BINDING, velocityBuffer);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, PARTICLE_FORCES_SSBO_BINDING, forcesBuffer);
    //glBindBufferBase(GL_SHADER_STORAGE_BUFFER, PARTICLE_DATA_SSBO_BINDING, particleDataBuffer);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, PARTICLE_SYSTEM_DATA_SSBO_BINDING, particleSystemDataBuffer);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, PARTICLE_EDGE_DATA_SSBO_BINDING, edgesBuffer);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, PARTICLE_HINGE_DATA_SSBO_BINDING, hingesBuffer);
    glBindBufferBase(GL_UNIFORM_BUFFER, PICKED_DATA_UBO_BINDING, pickedDebugDataBuffer);
    particleComputeShader->use();
}

void ClothSystem::SetupRender() {

    particleShader->use();
}

void ClothSystem::RenderDebug(GLuint VAO) {
    glBindVertexArray(VAO);
    glNamedBufferSubData(pickedDebugDataBuffer, 0, sizeof(PickedClothData), &pickedClothData);

    edgeDebugShader->use();
    glVertexArrayVertexBuffer(VAO, 0, edgeDebugDataBuffer, 0, sizeof(glm::vec4));
    glVertexArrayElementBuffer(VAO, edgeDebugEBO);
    glPointSize(10.0f);
    glLineWidth(10.0f);
    glDrawElements(GL_LINES, 2, GL_UNSIGNED_INT, 0);
    glDrawElements(GL_POINTS, 2, GL_UNSIGNED_INT, 0);

    hingeDebugShader->use();
    glPointSize(20.0f);
    glLineWidth(20.0f);
    glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
    glDrawElements(GL_POINTS, 6, GL_UNSIGNED_INT, 0);

    auto currTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<float, std::milli> duration = currTime - lastTime;
    
    if (duration.count() > 120) {
        pickedClothData.pickedEdgeIdx = (pickedClothData.pickedEdgeIdx + 1) % edgeCount;
        // fmt::print("Edge {}: {} {}\n",
        //     pickedClothData.pickedEdgeIdx, 
        //     edgeVec[pickedClothData.pickedEdgeIdx].x, 
        //     edgeVec[pickedClothData.pickedEdgeIdx].y);

        // fmt::print("Hinge {}: {} {} {} {}\n",
        //     pickedClothData.pickedHingeIdx, 
        //     hingeVec[pickedClothData.pickedHingeIdx].x, 
        //     hingeVec[pickedClothData.pickedHingeIdx].y,
        //     hingeVec[pickedClothData.pickedHingeIdx].z,
        //     hingeVec[pickedClothData.pickedHingeIdx].w);
        pickedClothData.pickedHingeIdx = (pickedClothData.pickedHingeIdx + 1) % hingeCount;        
        lastTime = currTime;
    }


}