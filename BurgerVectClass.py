# Class uses GTensor to evaluate Burgers Vector from DAXM DataSet

#!/usr/bin/env python 

import numpy as np
import os
import sys
from numpy import linalg as La

class BurgerVectClass:
  """Initialize Arrays corresponding the the Burger Vector and Slip Plane Normal Variants.
  """

  def __init__(self,BV,nMB,c_a):
    """ Initialize """
    self.BV=BV  # Initialize array for Burgers vectors 
    self.nMB=nMB # Initialize arrays for Plane Normals
    self.c_a=c_a # c to a ratio
    
  def BVCar(self):
    """ Convert Burgers Vectors to Cartesian Equivalent Coordinates. """
    BVec=self.BV
    BVCart=np.zeros((len(BVec),3))
    for ii in range(len(BVec)):
      BVCart[ii,0]=BVec[ii,0]
      BVCart[ii,1]=(BVec[ii,0]+2*BVec[ii,1])/np.sqrt(3.0)
      BVCart[ii,2]=BVec[ii,-1]/self.c_a
      BVCart[ii,:]=BVCart[ii,:]/La.norm(BVCart[ii,:])
    return BVCart
    
  def PlaneCart(self):
    """ Convert Slip Plane Normals to Cartesian Equivalets. """
    nplanes=len(self.nMB)
    PlaneNormal=np.zeros((nplanes,3)) 

    for ii in range(nplanes):
      PlaneNormal[ii,0]=self.nMB[ii,0]
      PlaneNormal[ii,1]=(self.nMB[ii,0]+2*self.nMB[ii,1])/np.sqrt(3.0)
      PlaneNormal[ii,2]=self.nMB[ii,-1]/self.c_a
      PlaneNormal[ii,:]=PlaneNormal[ii,:]/La.norm(PlaneNormal[ii,:])
    return PlaneNormal

    
  def GTn(self, GT,PlaneNormal):
    """ Perform GTn Operation to determine Burgers Vectors of dislocations piercing a plane pi with plane normal n """
    Index=[]
    npln=len(PlaneNormal)
    BurgerVec=[]
    for ii in range(len(GT)):
      if La.norm(GT[ii])!=0:
        GT[ii,:,:]=GT[ii,:,:]/La.norm(GT[ii,:,:])
        Index.append(ii)
        for jj in range(npln):
          BurgerVec.append(np.dot(GT[ii,:,:].T,PlaneNormal[jj,:]))
    BurgerVec=np.asarray([BurgerVec[ii]/La.norm(BurgerVec[ii]) for ii in range(len(BurgerVec))])
    return BurgerVec, npln*Index

 
  def BurgerVec(self,BVList,BVCart,ID,Tol):
    """ Method to Isolate Burgers Vectors from GTn based on Comparison with Standard List of Burgers Vectors."""
    Tol=np.cos(np.pi*5/180.) # 5-degree Tolerance for Vectors
    Bin,Index=[],[]
    for ii in range(len(BVList)):
      if np.abs(np.dot(BVCart,BVList[ii]))>=Tol:
        Bin.append(BVList[ii,:])
        Index.append(ID[ii])
    #Index=np.asarray(Index)
    return Bin,Index   # returns a list of the Burgers vectors that match BVCart and a non-unique list of the point indices

