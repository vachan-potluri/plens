# assumes initial diapgragm at x=0
import numpy as np

def shockTubeSoln(pLeft, pRight, TLeft, TRight, uLeft, uRight, xVec, time, gamma=1.4, R=287.0):
    xVec=np.asarray(xVec)
    # Following Knights book
    def pressureFn(p, pRef):
        # returns f(p, pRef). See eqn (2.98)
        if p>pRef:
            return 1/gamma*(p/pRef-1)/np.sqrt( 1/(2*gamma) * ( (gamma+1)*p/pRef + gamma-1 ) )
        else:
            return 2/(gamma-1)*( pow( p/pRef, (0.5*(gamma-1)/gamma) )-1 )
    
    def pressureFnDiff(p, pRef):
        # returns derivative of f(p, pRef) w.r.t p. See eqn (2.98)
        if p>pRef:
            return 1/pRef*np.sqrt(2/gamma)*( 1/np.sqrt((gamma+1)*p/pRef+gamma-1) -\
                                             0.5*(p/pRef-1)*pow((gamma+1)*p/pRef+gamma-1,-1.5)*(gamma+1))
        else:
            return 1/gamma*pow(p/pRef, -(0.5+0.5/gamma))*1/pRef
    
    aLeft=np.sqrt(gamma*R*TLeft)
    aRight=np.sqrt(gamma*R*TRight)
    rhoLeft=pLeft/(R*TLeft)
    rhoRight=pRight/(R*TRight)
    
    # initial guess for pressure in star region
    pStarNew=(rhoLeft*aLeft*pRight + rhoRight*aRight*pLeft + rhoLeft*rhoRight*aLeft*aRight*(uLeft-uRight))/(rhoLeft*aLeft + rhoRight*aRight)
    if pStarNew <= 1e-4:
        # very low pressure
        pStarNew = 1e-4
    pStarOld=-10000
    
    # solving for pStar
    print("Solving for pStar")
    while(abs(pStarOld-pStarNew) >= 0.001*pStarOld):
        pStarOld = pStarNew
        numer = aLeft*pressureFn(pStarOld,pLeft) + aRight*pressureFn(pStarOld,pRight) -uLeft + uRight
        denom = aLeft*pressureFnDiff(pStarOld,pLeft) + aRight*pressureFnDiff(pStarOld,pRight)
        pStarNew = pStarOld - numer/denom
        print("pStar residual:",abs(pStarNew-pStarOld)/pStarNew*100,"%")
    pStar=pStarNew
    print("pStar =",pStar)

    # The four cases
    if(pStar > pLeft and pStar > pRight):
        leftSpeed = uLeft - aLeft*np.sqrt( (0.5+0.5/gamma)*(pStar/pLeft-1)+1 )
        rightSpeed = uRight + aRight*np.sqrt( (0.5+0.5/gamma)*(pStar/pRight-1)+1 )
        contactSpeed = uLeft - aLeft/gamma*(pStar/pLeft-1)/np.sqrt( (0.5+0.5/gamma)*pStar/pLeft + 0.5-0.5/gamma )
        
        print("Case 1")
        i=0
        pVec=xVec.copy()
        TVec=xVec.copy()
        uVec=xVec.copy()
        while(i < xVec.size):
            x=xVec[i]
            if (x < leftSpeed*time):
                pVec[i]=pLeft
                TVec[i]=TLeft
                uVec[i]=uLeft
            elif (x >= leftSpeed*time and x < contactSpeed*time):
                pVec[i]=pStar
                uVec[i]=contactSpeed
                TVec[i]=pStar/( R*rhoLeft*(gamma-1+(gamma+1)*pStar/pLeft)/(gamma+1+(gamma-1)*pStar/pLeft) )
            elif (x >= contactSpeed*time and x < rightSpeed*time):
                pVec[i]=pStar
                uVec[i]=contactSpeed
                TVec[i]=pStar/( R*rhoRight*(gamma-1+(gamma+1)*pStar/pRight)/(gamma+1+(gamma-1)*pStar/pRight) )
            else:
                pVec[i]=pRight
                TVec[i]=TRight
                uVec[i]=uRight
            i=i+1
    elif(pStar > pLeft and pStar <= pRight):
        # left moving shock, right moving expansion
        shockSpeed = uLeft - aLeft*np.sqrt( (0.5+0.5/gamma)*(pStar/pLeft-1)+1 )
        expHeadSpeed = uRight + aRight
        expTailSpeed = uRight - aRight*( 2/(gamma-1)-(gamma+1)/(gamma-1)*pow(pStar/pRight, 0.5-0.5/gamma) )
        contactSpeed = uLeft - aLeft/gamma*(pStar/pLeft-1)/np.sqrt( (0.5+0.5/gamma)*pStar/pLeft + 0.5-0.5/gamma )
        
        print("Case 2")
        i=0
        pVec=xVec.copy()
        TVec=xVec.copy()
        uVec=xVec.copy()
        while(i < xVec.size):
            x=xVec[i]
            if(x < shockSpeed*time):
                pVec[i]=pLeft
                TVec[i]=TLeft
                uVec[i]=uLeft
            elif(x >= shockSpeed*time and x < contactSpeed*time):
                pVec[i]=pStar
                uVec[i]=contactSpeed
                TVec[i]=pStar/( R*rhoLeft*(gamma-1+(gamma+1)*pStar/pLeft)/(gamma+1+(gamma-1)*pStar/pLeft) )
            elif(x >= contactSpeed*time and x < expTailSpeed*time):
                pVec[i]=pStar
                uVec[i]=contactSpeed
                TVec[i]=pStar/( R*rhoRight*pow(pStar/pRight, 1/gamma) )
            elif(x >= expTailSpeed*time and x < expHeadSpeed*time):
                pVec[i]=pRight*( pow( (gamma-1)/(gamma+1)/aRight*(x/time-uRight)+2/(gamma+1), 2*gamma/(gamma-1) ) )
                TVec[i]=pVec[i]/(R*rhoRight*pow(pVec[i]/pRight, 1/gamma))
                uVec[i]=2/(gamma+1)*( x/time-aRight+0.5*(gamma-1)*uRight )
            else:
                pVec[i]=pRight
                TVec[i]=TRight
                uVec[i]=uRight
            i=i+1
    elif(pStar <= pLeft and pStar > pRight):
        # right moving shock, left moving expansion
        shockSpeed = uRight + aRight*np.sqrt( (0.5+0.5/gamma)*(pStar/pRight-1)+1 )
        expHeadSpeed = uLeft - aLeft
        expTailSpeed = uLeft + aLeft*( 2/(gamma-1)-(gamma+1)/(gamma-1)*pow(pStar/pLeft, 0.5-0.5/gamma) )
        contactSpeed = uRight + aRight/gamma*(pStar/pRight-1)/np.sqrt( (0.5+0.5/gamma)*pStar/pRight + 0.5-0.5/gamma )
        print("Case 3")
        
        i=0
        pVec=xVec.copy()
        TVec=xVec.copy()
        uVec=xVec.copy()
        while(i < xVec.size):
            x=xVec[i]
            if (x >= shockSpeed*time):
                pVec[i]=pRight
                TVec[i]=TRight
                uVec[i]=uRight
            elif (x < shockSpeed*time and x >= contactSpeed*time):
                pVec[i]=pStar
                TVec[i]=pStar/( R*rhoRight*(gamma-1+(gamma+1)*pStar/pRight)/(gamma+1+(gamma-1)*pStar/pRight) )
                uVec[i]=contactSpeed
            elif (x < contactSpeed*time and x >= expTailSpeed*time):
                pVec[i]=pStar
                TVec[i]=pStar/( R*rhoLeft*pow(pStar/pLeft, 1/gamma) )
                uVec[i]=contactSpeed
            elif (x < expTailSpeed*time and x >= expHeadSpeed*time):
                pVec[i]=pLeft*( pow( (gamma-1)/(gamma+1)/aLeft*(uLeft-x/time)+2/(gamma+1), 2*gamma/(gamma-1) ) )
                TVec[i]=pVec[i]/(R*rhoLeft*pow(pVec[i]/pLeft, 1/gamma))
                uVec[i]=2/(gamma+1)*( x/time+aLeft+0.5*(gamma-1)*uLeft )
            else:
                pVec[i]=pLeft
                TVec[i]=TLeft
                uVec[i]=uLeft
            i=i+1
    elif(pStar <= pLeft and pStar <= pRight):
        leftHead = uLeft - aLeft
        leftTail = uLeft + aLeft*( 2/(gamma-1)-(gamma+1)/(gamma-1)*pow(pStar/pLeft, 0.5-0.5/gamma) )
        rightHead = uRight + aRight
        rightTail = uRight - aRight*( 2/(gamma-1)-(gamma+1)/(gamma-1)*pow(pStar/pRight, 0.5-0.5/gamma) )
        contactSpeed = uLeft + 2*aLeft/(gamma-1)*(1-pow(pStar/pLeft, 0.5-0.5/gamma))
        
        print("Case 4")
        i=0
        pVec=xVec.copy()
        TVec=xVec.copy()
        uVec=xVec.copy()
        while(i < xVec.size):
            x=xVec[i]
            if(x < leftHead*time):
                pVec[i]=pLeft
                TVec[i]=TLeft
                uVec[i]=uLeft
            elif(x >= leftHead*time and x < leftTail*time):
                pVec[i]=pLeft*( pow( (gamma-1)/(gamma+1)/aLeft*(uLeft-x/time)+2/(gamma+1), 2*gamma/(gamma-1) ) )
                TVec[i]=pVec[i]/(R*rhoLeft*pow(pVec[i]/pLeft, 1/gamma))
                uVec[i]=2/(gamma+1)*( x/time+aLeft+0.5*(gamma-1)*uLeft )
            elif(x >= leftTail*time and x < contactSpeed*time):
                pVec[i]=pStar
                TVec[i]=pStar/( R*rhoLeft*pow(pStar/pLeft, 1/gamma) )
                uVec[i]=contactSpeed
            elif(x >= contactSpeed*time and x < rightTail*time):
                pVec[i]=pStar
                uVec[i]=contactSpeed
                TVec[i]=pStar/( R*rhoRight*pow(pStar/pRight, 1/gamma) )
            elif(x >= rightTail*time and x < rightHead*time):
                pVec[i]=pRight*( pow( (gamma-1)/(gamma+1)/aRight*(x/time-uRight)+2/(gamma+1), 2*gamma/(gamma-1) ) )
                TVec[i]=pVec[i]/(R*rhoRight*pow(pVec[i]/pRight, 1/gamma))
                uVec[i]=2/(gamma+1)*( x/time-aRight+0.5*(gamma-1)*uRight )
            else:
                pVec[i]=pRight
                TVec[i]=TRight
                uVec[i]=uRight
            i=i+1
    return np.column_stack((xVec,TVec,pVec,uVec,pVec/(R*TVec)))
