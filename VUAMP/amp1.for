!     user amplitude subroutine
      Subroutine VUAMP(
!          passed in for information and state variables
     *     ampName, time, ampValueOld, dt, nprops, props, nSvars, 
     *     svars, lFlagsInfo, nSensor, sensorValues, sensorNames, 
     *     jSensorLookUpTable, 
!          to be defined
     *     ampValueNew, 
     *     lFlagsDefine,
     *     AmpDerivative, AmpSecDerivative, AmpIncIntegral)
      
      include 'vaba_param.inc'

!     svars - additional state variables, similar to (V)UEL
      dimension sensorValues(nSensor), props(nprops), 
     *     svars(nSvars)
      character*80 sensorNames(nSensor)
      character*80 ampName

!     time indices
      parameter( iStepTime        = 1,
     *           iTotalTime       = 2,
     *           nTime            = 2)
!     flags passed in for information
      parameter( iInitialization   = 1,
     *           iRegularInc       = 2,
     *           ikStep            = 3,
     *           nFlagsInfo        = 3)
!     optional flags to be defined
      parameter( iComputeDeriv     = 1,
     *           iComputeSecDeriv  = 2,
     *           iComputeInteg     = 3,
     *           iStopAnalysis     = 4,
     *           iConcludeStep     = 5,
     *           nFlagsDefine      = 5)

      parameter( tStep=0.18, tAccelerateMotor = .00375, 
     *           omegaFinal=23.26)

!				Alternatively, assign the user-defined amplitude 
!				properties on the data lines rather than using a parameter
!				definition above.
!				 tStep = props(1)
!				 tAccelerateMotor = props(2)
!				 omegaFinal = props(3)      

      dimension time(nTime), lFlagsInfo(nFlagsInfo),
     *          lFlagsDefine(nFlagsDefine)
      dimension jSensorLookUpTable(*)

      lFlagsDefine(iComputeDeriv)    = 1
      lFlagsDefine(iComputeSecDeriv) = 1

!     get sensor value
      vTrans_CU1  = vGetSensorValue('HORIZ_TRANSL_MOTION',
     *                              jSensorLookUpTable,
     *                              sensorValues)

      if (ampName(1:22) .eq. 'MOTOR_WITH_STOP_SENSOR' ) then
         if (lFlagsInfo(iInitialization).eq.1) then 
            ampValueNew      = ampValueOld

            svars(1) = 0.0
            svars(2) = 0.0
         else
            tim = time(iStepTime)

!           ramp up the angular rot velocity  of the electri! 
!           motor after which hold constant
            if (tim .le. tAccelerateMotor) then 
               ampValueNew = omegaFinal*tim/tAccelerateMotor

            else
               ampValueNew = omegaFinal
            end if

!           retrieve old sensor value
            vTrans_CU1_old = svars(1)

!           detect a zero crossing and count the number of 
!           crossings
            if (vTrans_CU1_old*vTrans_CU1 .le. 0.0 .and.
     *          tim .gt. tAccelerateMotor ) then 
               svars(2) = svars(2) + 1.0
            end if            
            nrCrossings = int(svars(2))

!           stop the motor if sensor crosses zero the second 
!           time
            if (nrCrossings.eq.2) then 
               ampValueNew =  0.0
               lFlagsDefine(iConcludeStep)=1
            end if

!           store sensor value
            svars(1) = vTrans_CU1

         end if         
      end if 

      return
      end