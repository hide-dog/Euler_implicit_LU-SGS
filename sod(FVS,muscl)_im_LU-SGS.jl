using Printf

# --------------------------
# -- function setup       --
# --------------------------
#
# 基本量と保存量の設定
#
function setup(nx,dx,gamma)
    u   = zeros(nx)               # 速度
    rho = zeros(nx)               # 密度
    p   = zeros(nx)               # 圧力
    e   = zeros(nx)               # エネルギー
    x   = zeros(nx)               # 位置

    qf = zeros(nx,3)              # 基本量
    Qc = zeros(nx,3)              # 保存量


    """
    e=p/(r-1)+rho*u^2/2
    """
    for i in 1:nx
        u[i] = 0.0

        if i <= nx*0.5
            rho[i] = 1.0
            p[i] = 1.0
        elseif i >= nx*0.5
            rho[i] = 0.125
            p[i] = 0.1
        end
        e[i] = p[i]/(gamma-1)+rho[i]*(u[i]^2)/2
        x[i] = i*dx-dx/2
    end


    for i in 1:nx
        qf[i,1] = u[i]
        Qc[i,1] = rho[i]
    
        qf[i,2] = rho[i]
        Qc[i,2] = u[i]*rho[i]
    
        qf[i,3] = p[i]
        Qc[i,3] = e[i]
    
    end
    return x,qf,Qc
end

# --------------------------
# -- function boundary    --
# --------------------------
#
# 境界条件の設定
#
function boundary(nx,gamma,qf,Qc,bdcon)
    if Int(bdcon[1][1]) == 0
        for i in 1:3
            Qc[nx,i] = Qc[nx-1,i]
        end
    elseif Int(bdcon[1][1]) == 1
        rho = bdcon[1][2]
        u = bdcon[1][3]
        p = bdcon[1][4]
        e = p/(gamma-1)+rho*(u^2)/2
        Qc[nx,1] = rho
        Qc[nx,2] = rho*u
        Qc[nx,3] = e
    end

    if Int(bdcon[1][2]) == 0
        for i in 1:3
            Qc[1,i] = Qc[2,i]
        end
    elseif Int(bdcon[1][2]) == 1
        rho = bdcon[2][2]
        u = bdcon[2][3]
        p = bdcon[2][4]
        e = p/(gamma-1)+rho*(u^2)/2
        Qc[1,1] = rho
        Qc[1,2] = rho*u
        Qc[1,3] = e
    end

    return Qc
end
        
# ---------------------------
# -- function LDU_forlusgs --
# ---------------------------
#
# LU-SGSで行うLDU分解に向けて、ヤコビアン行列の近似とbetaの計算を行う
#
function LDU_forlusgs(Qc,qf,nx,beta,gamma)
    Amatrix_plus  = zeros(nx,3,3)
    Amatrix_minus = zeros(nx,3,3)
    beta_sigma    = zeros(nx)

    I = zeros(3,3)
    for i in 1:3
        I[i,i] = 1.0
    end

    for i in 1:nx
        H = (Qc[i,3]+qf[i,3])/Qc[i,1]  # エンタルピー
        u = qf[i,1]
        c = ((gamma-1)*(H-0.5*u^2))^0.5
        b_para = (gamma-1)/c^2
        a_para = 0.5*b_para*u^2

        sigma = abs(u) + c
        # iセルにおけるR,R^-1,Λ,|Λ|
        R, R_inv, Gam, Gam_abs = A_pm(H,u,c,b_para,a_para)
        A_matrix = nn_inner_product((nn_inner_product(R, Gam)), R_inv)
        
        for j in 1:3
            for k in 1:3
                temp = I[j,k]*beta*sigma
                Amatrix_plus[i,j,k] = 0.5*(A_matrix[j,k]+temp)
                Amatrix_minus[i,j,k] = 0.5*(A_matrix[j,k]-temp)
            end
        end
        beta_sigma[i] = beta*sigma
    end
    
    return Amatrix_plus,Amatrix_minus,beta_sigma
end

# --------------------------
# -- function cal_RHS     --
# --------------------------
#
# 右辺の計算
# RHSを定義
#
function cal_RHS(Fplus,nx)  # 境界フラックスの計算

    RHS = zeros(nx,3)

    for i in 2:nx-1
        for j in 1:3
            RHS[i,j] = Fplus[i,j]-Fplus[i-1,j]
        end
    end
    return RHS
end

# --------------------------
# -- function fvs         --
# --------------------------
#
# Flux Vector Splitting method 
# Fplusを定義
#
function fvs(QcL,QcR,qfL,qfR,nx,gamma)  # FVS法によるフラックスの計算(セル1と2の境界をFplus[1]に格納)

    Fplus = zeros(nx+1,3)

    for i in 1:nx-1
        H = (QcL[i,3]+qfL[i,3])/QcL[i,1]  # エンタルピー
        u = qfL[i,1]
        c = ((gamma-1)*(H-0.5*u^2))^0.5
        b_para = (gamma-1)/c^2
        a_para = 0.5*b_para*u^2
        # i+1/2セルにおけるR,R^-1,Λ,|Λ|
        R, R_inv, Gam, Gam_abs = A_pm(H,u,c,b_para,a_para)
        Ap = nn_inner_product((nn_inner_product(R, Gam+Gam_abs)), R_inv)               # 固有値が正のものを計算
        
        H = (QcR[i,3]+qfR[i,3])/QcR[i,1]  # エンタルピー
        u = qfR[i,1]
        c = ((gamma-1)*(H-0.5*u^2))^0.5
        b_para = (gamma-1)/c^2
        a_para = 0.5*b_para*u^2
        # i+1/2セルにおけるR,R^-1,Λ,|Λ|
        R, R_inv, Gam, Gam_abs = A_pm(H,u,c,b_para,a_para)
        Am = nn_inner_product((nn_inner_product(R, Gam-Gam_abs)), R_inv)               # 固有値が負のものを計算
        
        tempL = zeros(3)
        tempR = zeros(3)
        for j in 1:3
            tempL[j] = QcL[i,j]
            tempR[j] = QcR[i,j]
        end

        tempL = nm_inner_product(Ap, tempL)
        tempR = nm_inner_product(Am, tempR)

        for j in 1:3
            Fplus[i,j] = 0.5*(tempL[j] + tempR[j])   # フラックスを計算tempL[j] = 
        end
    end
    return Fplus
end

# --------------------------
# -- function A_pm       --
# --------------------------
#
# ヤコビアン行列の正負の計算に向け、固有値等を計算
# A = R^(-1)ΛR としたとき
# R = R
# R_inv = R^(-1)
# Gam = Λ
#
function A_pm(H,u,c,b_para,a_para)  # ヤコビアン行列の固有値もろもろ計算
    R     = [1.0 1.0 1.0;
            u-c u u+c;
            H-u*c 0.5*u^2 H+u*c]
    R_inv = [0.5*(a_para+u/c) 0.5*(-b_para*u-1/c) 0.5*b_para;
            1-a_para b_para*u -b_para;
            0.5*(a_para-u/c) 0.5*(-b_para*u+1/c) 0.5*b_para]
    Gam   = [(u-c) 0.0 0.0; 0.0 u 0.0; 0.0 0.0 (u+c)]
    Gam_abs = [abs(u-c) 0.0 0.0; 0.0 abs(u) 0.0; 0.0 0.0 abs(u+c)]

    return R, R_inv, Gam, Gam_abs
end

# --------------------------
# -- function yacobi_A    --
# --------------------------
#
# ヤコビアン行列の計算
#
function yacobi_A(Qc,qf,gamma,nx)
    yacobiA = zeros(nx,4,4)
    
    for i in 1:nx
        H = (Qc[i,3]+qf[i,3])/Qc[i,1]  # エンタルピー
        u = qf[i,1]
        c = ((gamma-1)*(H-0.5*u^2))^0.5
        b_para = (gamma-1)/c^2
        a_para = 0.5*b_para*u^2

        R, R_inv, Gam, Gam_abs = A_pm(H,u,c,b_para,a_para)

        temp = nn_inner_product((nn_inner_product(R, Gam)), R_inv)

        for j in 1:3
            for k in 1:3
                yacobiA[i,j,k] = temp[j,k]
            end
        end
    end

    return yacobiA
end

# -------------------------------
# -- function nn_inner_product --
# -------------------------------
#
# n*n行列同士の内積
#
function nn_inner_product(a,b)
    temp=zeros(size(a)[1],size(a)[1])
    for i in 1:size(a)[1]
        for j in 1:size(a)[1]
            for k in 1:size(a)[1]
                temp[i,j] += a[i,k]*b[k,j]
            end
        end
    end
    return temp    
end

# -------------------------------
# -- function nm_inner_product --
# -------------------------------
#
# n*m行列と1*m行列の内積
#
function nm_inner_product(a,b)
    temp = zeros(size(a)[1])
    for i in 1:size(a)[1]
        for k in 1:size(a)[2]
            temp[i] += a[i,k]*b[k]
        end
    end
    return temp
end

# -------------------------------
# -- function muscl            --
# -------------------------------
#
# muscl法による補間
#
function muscl(Qc,qf,nx,k_muscl,b_muscl,gamma)

    qfL = zeros(nx+1,3)
    qfR = zeros(nx+1,3)

    for i in 2:nx-2
        for j in 1:3
            dplus_j   = qf[i+1,j]-qf[i,j]
            dminus_j  = qf[i,j]-qf[i-1,j]
            dplus_jp  = qf[i+2,j]-qf[i+1,j]
            dminus_jp = qf[i+1,j]-qf[i,j]
            
            qfL[i,j] = qf[i,j]+1/4*((1-k_muscl)*minmod(dminus_j,dplus_j,b_muscl)+(1+k_muscl)*minmod(dplus_j,dminus_j,b_muscl))
            qfR[i,j] = qf[i+1,j]-1/4*((1-k_muscl)*minmod(dplus_jp,dminus_jp,b_muscl)+(1+k_muscl)*minmod(dminus_jp,dplus_jp,b_muscl))
        end
    end

    # 境界内側用
    for j in 1:3
        dplus_jp  = qf[3,j]-qf[2,j]
        dminus_jp = qf[2,j]-qf[1,j]
        qfR[1,j] = qf[2,j]-1/4*((1-k_muscl)*minmod(dplus_jp,dminus_jp,b_muscl)+(1+k_muscl)*minmod(dminus_jp,dplus_jp,b_muscl))

        dplus_j      = qf[nx,j]-qf[nx-1,j]
        dminus_j     = qf[nx-1,j]-qf[nx-2,j]
        qfL[nx-1,j] = qf[nx-1,j]+1/4*((1-k_muscl)*minmod(dminus_j,dplus_j,b_muscl)+(1+k_muscl)*minmod(dplus_j,dminus_j,b_muscl))
    end
        
            
    QcL = qftoQc(qfL,nx,gamma)
    QcR = qftoQc(qfR,nx,gamma)

    # 境界外側用(境界は風上)
    for i in 1:3
        qfL[1,i] = qf[1,i]
        QcL[1,i] = Qc[1,i]
        qfR[nx-1,i] = qf[nx,i]
        QcR[nx-1,i] = Qc[nx,i]
    end

    return QcL,QcR,qfL,qfR
end
    
# -------------------------------
# -- function minmod           --
# -------------------------------
#
# 流速制限関数minmod
#
function minmod(x,y,b)
    ans = sign(x)*max(0,min(abs(x),sign(x)*y*b))
    return ans
end

# -------------------------------
# -- function qftoQc           --
# -------------------------------
#
# 基本量から保存量に変換
# 
function qftoQc(qf,nx,gamma)  # 基本量から保存量変換
    lo_Qc = zeros(nx,3)
    for i in 1:nx
        lo_Qc[i,1] = qf[i,2]
        lo_Qc[i,2] = qf[i,2]*qf[i,1]
        lo_Qc[i,3] = (qf[i,3]/(gamma-1)+1.0/2.0*qf[i,2]*(qf[i,1]^2))
    end

    return lo_Qc
end

# -------------------------------
# -- function Qctoqf           --
# -------------------------------
#
# 保存量から基本量に変換
#
function Qctoqf(Qc,nx,gamma)  # 保存量から基本量変換
    lo_qf = zeros(nx,3)
    for i in 1:nx
        lo_qf[i,1] = Qc[i,2]/Qc[i,1]
        lo_qf[i,2] = Qc[i,1]
        lo_qf[i,3] = (gamma-1)*(Qc[i,3]-1.0/2.0*Qc[i,1]*((Qc[i,2]/Qc[i,1])^2))
    end
    return lo_qf
end

# -------------------------------
# -- function lusgs        --
# -------------------------------
#
# LU-SGSによる反復
#
function lusgs(Qcn,Qcm,nx,dx,dt,RHS,Amatrix_plus,Amatrix_minus,beta_sigma,norm_ok) # 内部反復
    delta_Q = zeros(nx,3)
    delta_Q_temp = zeros(nx,3)

    sum_b = zeros(3)
    for i in 2:nx-1
        for j in 1:3
            sum_b[j] += abs(RHS[i,j])
        end
    end

    ite=0
    con=0
    while con == 0 # lusgs loop
        delta_Q_temp = copy(delta_Q)

        # LDU分解
        L = zeros(nx,3)
        D = zeros(nx)
        U = zeros(nx,3)

        for i in 1:nx
            temp_plus = zeros(3,3)
            temp_minus = zeros(3,3)
            temp = zeros(3)
            for l in 1:3
                for m in 1:3
                    temp_plus[l,m] = Amatrix_plus[i,l,m]
                    temp_minus[l,m] = Amatrix_minus[i,l,m]
                end
                temp[l] = delta_Q[i,l]
            end
            temp_plus = nm_inner_product(temp_plus,temp)
            temp_minus = nm_inner_product(temp_minus,temp)

            for j in 1:3
                L[i,j] = dt*0.5*temp_plus[j]
                U[i,j] = dt*0.5*temp_minus[j]
            end
            D[i] = dx+dx+dt*(beta_sigma[i])
        end

        RHS_temp = zeros(nx,3)
        # RHS
        for i in 2:nx-1
            for j in 1:3
                RHS_temp[i,j] = -(Qcm[i,j]-Qcn[i,j])*dx-dt*RHS[i,j]
            end
        end

        # (D+L)Q=Q
        for i in 2:nx-1
            for j in 1:3
                delta_Q[i,j] = (L[i-1,j]+RHS_temp[i,j]) / D[i]
            end
        end

        # (D+U)Q=DQ
        for i in 2:nx-1
            for j in 1:3
                delta_Q[i,j] = delta_Q[i,j] - U[i+1,j]/ D[i]
            end
        end
        
        # 収束判定
        if (ite+1) % 100 == 0
            sum_b_Ax = zeros(3)

            for i in 2:nx-1
                for j in 1:3
                    sum_b_Ax[j] += abs(delta_Q[i,j]-delta_Q_temp[i,j])
                end
            end
            
            norm2d = zeros(3)

            for i in 1:3
                norm2d[i] = sum_b_Ax[i]/sum_b[i]
            end

            if norm2d[1] < norm_ok && norm2d[2] < norm_ok && norm2d[3] < norm_ok
                con=1
            end
        end
        ite += 1

        if ite > 1.0e5
            println("Number of iterations exceeded 1.0e5 ")
            throw(UndefVarError(:x))
        end
    end
    
    # Qcの更新
    for i in 2:nx-1
        for j in 1:3
            Qcm[i,j] = Qcm[i,j]+delta_Q[i,j]
        end
    end
    return Qcm
end

# -------------------------------
# -- function output_q         --
# -------------------------------
#
# x,rho,u,pの出力
#
function output_q(x,qf,nx,out_dir,out_file_front,stepnum,out_ext)
    fff = out_dir*"/"*out_file_front*string(stepnum)*out_ext
    open(fff,"w") do f
        write(f,"result:x, rho[kg/m^3], u[m/s], p[Pa]\n")
        for i in 2:nx-1
            a = @sprintf("%8.8e", x[i])
            b = @sprintf("%8.8e", qf[i,1])
            c = @sprintf("%8.8e", qf[i,2])
            d = @sprintf("%8.8e", qf[i,3])
            write(f,a*" "*b*" "*c*" "*d*"\n")
        end
    end
end

# -------------------------------
# -- function cre_dir          --
# -------------------------------
#
# ディレクトリの作成
#
function cre_dir(outdir)
    k = 0
    try rm(outdir,recursive=true)
    catch
        mkdir(outdir)
        k = 1
    end

    if k == 0
        mkdir(outdir)
    end
end


# -------------------------------
# -- main                      --
# -------------------------------
function main()
    # --------------------------
    # -- initial value        --
    # --------------------------
    nstep = 300                     # 時間ステップ数
    nx0   = 100                     # 空間ステップ数
    dt    = 0.002                   # 時間刻み幅
    dx    = 0.01                    # 空間刻み幅
    innner_step = 10
    dtau = 0.002

    # -- slope limiter --
    k_muscl = 1/3                      # muscl精度
    b_muscl = (3-k_muscl)/(1-k_muscl)

    # -- constant --
    gamma = 1.4                    # 比熱比

    # -- boundary --
    # bd1_con = 0 :流入条件
    #         = 1 :流出条件
    #
    # -- x+ right --
    bd1_con = 0
    bd1_rho = 1.0
    bd1_u   = 0.0
    bd1_p   = 1.0
    bd1_T   = 300.0
    # -- x- left --
    bd2_con = 1
    bd2_rho = 1.0
    bd2_u   = 0.0
    bd2_p   = 1.0
    bd2_T   = 300.0

    # -- itetarive --
    norm_ok = 1.0e-4               # 収束条件
    beta = 1.1

    # -- 出力--
    dir_name       = "sod_lusgs_0.002s_jl"           # 出力フォルダ名
    out_name_front = "time"                 # 出力ファイル名（先頭）
    out_name_back  = "d-3"

    # --------------------------
    # -- prepare cal          --
    # --------------------------
    lbound = 1                     # 仮想境界セル数
    nx = nx0+2*lbound              # 総空間セル数

    bdcon = [[bd1_con, bd1_rho, bd1_u, bd1_p, bd1_T],
            [bd2_con, bd2_rho, bd2_u, bd2_p, bd2_T]]
    
    # --------------------------
    # -- setup                --
    # --------------------------
    cre_dir(dir_name)
    x,qf,Qc = setup(nx,dx,gamma)

    # --------------------------
    # -- main loop            --
    # --------------------------
    for t in 1:nstep
        print(t)
        
        Qcn = copy(Qc)
        Qcm = copy(Qc)

        # 内部反復
        for innner_t in 1:innner_step
            
            qfm = Qctoqf(Qcm,nx,gamma)
            Qcm = boundary(nx,gamma,qfm,Qcm,bdcon)
            QcL,QcR,qfL,qfR = muscl(Qcm,qfm,nx,k_muscl,b_muscl,gamma)

            Fplus = fvs(QcL,QcR,qfL,qfR,nx,gamma) 
            RHS = cal_RHS(Fplus,nx)

            Amatrix_plus,Amatrix_minus,beta_sigma = LDU_forlusgs(Qcm,qfm,nx,beta,gamma)

            Qcm = lusgs(Qcn,Qcm,nx,dx,dt,RHS,Amatrix_plus,Amatrix_minus,beta_sigma,norm_ok)
        end

        Qc = copy(Qcm)    
        qf = Qctoqf(Qc,nx,gamma)

        stepnum = Int(round(t*dt*1000))
        output_q(x,qf,nx,dir_name,out_name_front,stepnum,out_name_back)

        #throw(UndefVarError(:x))
    end
end

# main
main()



