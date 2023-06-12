      Subroutine VUMAT(
     1  nblock,ndir,nshr,nstatev,nfieldv,nprops,lanneal,stepTime,
     2  totalTime,dt,cmname,coordMp,charLength,props,density,
     3  strainlnc,relSpinInc,tempOld,stretchOld,defgradOld,fieldold,
     4  stressOld,stateOld,enerInternOld,enerInelasOld,tempNew,
     5  stretchNew,defgradNew,fieldNew,dstrain
     6  stressNew,stateNew,enerInternNew,enerInelasNew)

        Implicit Double Precision (a-h, o-z)
        Parameter (j_sys_Dimension = 2, n_vec_Length = 136,
     1    maxblk = n_vec_Length)

        Double Precision :: props(nprops), density_abq(nblock), 
     1    coordMp(nblock,*), charLength(nblock,*),
     2    strainInc(nblock,ndir+nshr), relSpinInc(nblock,nshr), 
     3    tempOld(nblock), tempNew(nblock),
     4    stretchOld(nblock,ndir+nshr), 
     5    defgradOld(nblock,ndir+nshr+nshr), fieldOld(nblock,nfieldv),
     6    stressold(nblock,ndir+nshr), stateOld(nblock,nstatev),
     7    enerInternold(nblock), enerInelasOld(nblock),
     8    stretchNew(nblock,ndir+nshr), dstrain(nblock,ndir+nshr)
        
        Double Precision :: stressNew(nblock,ndir+nshr),
     1    defgradNew(nblock,ndir+nshr+nshr), fieldNew(nblock,nfieldv),
     2    stateNew(nblock,nstatev), enerInternNew(nblock),
     3    enerInelasNew(nblock)

        Double Precision :: stepTime, totalTime, dt

        Character(len=80) :: cmname

        Double Precision ep1, ep2, ep3, sb, sc, scc, k, st, stt, e1, e2,
     1    e3

        ! Read user—defined material properties
        ep1 = props(1) ! compressive strain at damage initiation
        ep2 = props(2) ! compressive strain after cell wall collapse
        ep3 = props(3) ! compressive strain after continued crushing
        sb = props(4) ! compressive strength
        sc = props(5) ! compressive stress after cell wall collapse
        scc = props(6) ! compressive stress after continued crushing
        k   = props(7) ! ratio of tensile modulus to compressive modulus
        st  = props(8) ! tensile strength
        stt = props(9) ! tensile stress at which a crushed core's tangent stiffness changes from compressive modulus

        ! Calculate dependent properties
        e1 = sb/ep1 ! compressive modulus
        e2 = (sb - sc)/(ep1 - ep2) ! tangent stiffness during cell wall collapse
        e3 = (scc - sc)/(ep3 - ep2) ! tangent stiffness during continued crushing

        ! --------------------------------------------------------------------- !
        !    Solution—dependent Variables (SDV):
        !    SDV1: Damage state, where:
        !              (—l)=Tensile loaded only (loading or unloading)
        !               (1)=Elastically compressed (O>strain>=ep1)
        !               (2)=Partially crushed (ep1>strain>=ep2)
        !               (3)=Crushed (ep2>strain)
        !    (inactive) (4)=Unloaded from crushed state with stress<stt
        !    (inactive) (5)=Unloaded from crushed state with stress>stt
        !    SDV2: Max applied -ve strain
        !    SDV3: Stress at max applied -ve strain
        !    SDV4: Strain when stressNew(1):stt
        !    SDV5: Strain when stress reaches stt during reloading
        !    SDV6: Indicator for stress exceeding stt, where:
        !               (O)=Stress didn't exceed stt
        !               (1)=Stress exceeded stt
        !    SDV7: Previous strain
        ! --------------------------------------------------------------------- !
        
        Master: Do km = 1,nblock

          Call CoreCrush(stressOld(km,1),stressNew(km,1),stateOld(km,:),
     1                    stateNew(km,:),dstrain(km,1),nstatev,ep1,
     2                    ep2,ep3,sb,sc,scc,k,st,stt,e1,e2,e3)
          
                          
        end Do Master
      ! --------------------------------------------------------------------- !  
        Return
      end Subroutine VUMAT

      
      Subroutine CoreCrush(stressOld,stressNew,stateOld,stateNew,
     1                    dstrain,nstatev,ep1,ep2,ep3,sb,sc,scc,k,
     2                    st,stt,e1,e2,e3)

        Integer, intent(IN) :: nstatev
        Double Precision, intent(IN) :: stressOld, stateOld(nstateV)
        Double Precision, intent(OUT) :: stressNew, stateNew(nstateV)
        Double Precision, intent(IN) :: dstrain
        Double Precision, intent(IN) :: ep1,ep2,ep3,sb,sc,scc,k,st,stt,
     1    e1,e2,e3

        Double Precision :: DDSDDE, strainold, strain
      ! --------------------------------------------------------------------- !
        stateNew(:) = stateOld(:) ! Initialize state Variables

        strainOld = stateOld(7) ! Previous strain
        strain = strainold + dstrain ! Current strain
        stateNew(7) = strain ! Store updated strain

        If (strain > 0 .AND. stateOld(2) >= ep1) Then ! Tensile loading
          DDSDDE = k*e1
          stressNew = stressold + DDSDDE*dstrain
          stateNew(1) = -1

        Else If (dstrain <= 0) Then ! Compressive loading
          If (strain <= stateOld(2)) Then ! If strain is at a new high,
            Loading: If (strain >= ep1) Then ! Elastically Compressed
              DDSDDE = e1
              stateNew(1) = 1
            Else If (strain >= ep2) Then Loading ! Partially crushed
              DDSDDE = e2
              stateNew(1) = 2
            Else Loading ! Fully crushed
              DDSDDE = e3
              stateNew(1) = 3
            end If Loading

            stressNew = stressold + DDSDDE*dstrain
            stateNew(2) = strain
            stateNew(3) = stressNew
            stateNew(6) = 0

          Else ! If the core is being reloaded,
            Reloading: If (stateOld(2) >= ep1) Then  ! Elastically compressed
              DDSDDE = e1
              stressNew = stressold + DDSDDE*dstrain
              stateNew(1) = 1

            Else If (stateOld(2) >= ep2) Then Reloading ! Partially Crushed
              DDSDDE = stateOld(3)/stateOld(2)
              stressNew = stressold + DDSDDE*dstrain

            Else Reloading ! Fully Crushed
              If (stateOld(6) == 1) Then ! Unloaded past stt
                If (stressOld >= stt) Then
                  DDSDDE = e1
                  stressNew = stressold + DDSDDE*dstrain
                  If (stressNew < stt) Then
                    stateNew(5) = (stt - stressold)/e1 + strainOld
                    DDSDDE = ((stt - stateOld(3))/(stateNew(5)
     1              - stateOld(2)))
                    stressNew = stt + DDSDDE*(strain - stateNew(5))
                  end If
                Else
                  DDSDDE = (stressOld - stateOld(3))/(strainold
     1            - stateOld(2))
                  stressNew = stressOld + DDSDDE*dstrain
                end If

              Else ! Not unloaded past stt
                If (stressOld < stt) Then
                  DDSDDE = e1
                  stressNew = stressOld + DDSDDE*dstrain
                end If
              end If
            end If Reloading
          end If

        Else ! If the core is being unloaded,
      !   Unloading when undamaged:
          Unloading: If (stateOld(2) >= ep1) Then
            DDSDDE = e1
            stressNew = stressold + DDSDDE*dstrain

      !   Unloding when partially Crushed:
          Else If (stateOld(2) >= ep2) Then Unloading
            DDSDDE = stateOld(3)/stateOld(2)
            stressNew = stressold + DDSDDE*dstrain
            If (stressNew >= stt) stateNew(6) = 1

      !   Unloading when fully crushed:
          Else Unloading
            If (stressOld <= stt .AND. stateOld(6) == 0) Then
              DDSDDE = e1
              stressNew = stressold + DDSDDE*dstrain
              If (stressNew > stt) Then
                stateNew(4) = strainOld + (stt - stressOld)/e1
                DDSDDE = -(st - stt)/stateNew(4)
                stressNew = stt + DDSDDE*(strain - stateNew(4))
                stateNew(6) = 1
              end If

            Else If (stressOld <= stt.AND.stateOld(6) == 1) Then
              DDSDDE = ((stt - stateOld(3))/(stateOld(5) - stateOld(2)))
              stressNew = stressold + DDSDDE*dstrain
              If (stressNew > stt) Then
                DDSDDE = e1
                stressNew = stt + DDSDDE*(strain - stateOld(5))
              end If

            Else
              DDSDDE = (st - stressold)/(-strainold)
              stressNew = stressold + DDSDDE*dstrain
              stateNew(6) = 1
            end If

          end If Unloading
        end If

        Return
      end subroutine Corecrush

