import numpy as np
import time

def TDMA(T_inf,T_b,k,D,L,h,alfa,delta_t,Tp_zero,n_volumes,criterio):
    A = ((np.pi)*D**2)/4
    P = (np.pi)*D
    alfa_k = alfa/k
    deltax = L/n_volumes
    #Mp = deltax
    A_e = 1/(deltax)
    A_w = 1/(deltax)
    deltax_fronteira = (deltax)/2
    A_fronteira = 1/(deltax_fronteira)
    Ap_esquerda = A_e + A_fronteira + (deltax*h*P/(A*k)) + deltax/(alfa*delta_t)
    Ap_meio = A_e + A_w + (deltax*h*P/(k*A)) + deltax/(alfa*delta_t)
    Ap_direita = A_w + deltax/(alfa*delta_t) + h*deltax*P/(k*A)
    
  

    Pm=np.zeros(n_volumes)
    Qm=np.zeros(n_volumes)
    Tm_i = np.zeros(n_volumes)
    Tm = np.zeros(n_volumes)
    B = np.zeros(n_volumes)
    
    # Start das temperaturas
    for i in range(0,n_volumes):
        Tm_i[i] = Tp_zero
        Tm[i] = Tp_zero
        
    erro = 1 # Start da código
    
    inicio = time.time()  #Calcular o tempo de início
    cont = 0

    
    while erro >criterio:
        
    #Cálculo do B
        for i in range(0,n_volumes):
            if i == 0:
                B[i] = (Tm_i[i]*deltax)/(alfa*delta_t) + (A_fronteira*T_b) + h*T_inf*deltax*P/(k*A)

            elif i == (n_volumes-1):
                B[i] = (h*T_inf*deltax*P/(k*A)) + (Tm_i[i]*deltax)/(alfa*delta_t)

            else:
                B[i] = (deltax*Tm_i[i])/(alfa*delta_t) + (h*T_inf*deltax*P/(k*A))


            #Calculo dos coeficientes do TDMA
            
            for i in range(0,n_volumes):
                #Calculo Esquerda
                if i == 0:
                    Pm[i] = A_e/Ap_esquerda
                    Qm[i] = B[i]/Ap_esquerda
                #Calculo Direita    
                elif i ==(n_volumes -1):
                    Qm[i] = (B[i]+A_w*Qm[i-1])/(Ap_direita-A_w*Pm[i-1])
                #Calculo meio   
                else:
                    Pm[i] = A_e/(Ap_meio-A_w*Pm[i-1])
                    Qm[i] = (B[i]+A_w*Qm[i-1])/(Ap_meio-A_w*Pm[i-1])


            #Calculo Temperatura
            for i in range(n_volumes-1,-1,-1):
                #Tempeatura direita
                if i == (n_volumes-1):
                    Tm[i] = Qm[i]

                #Outras Temperaturas    
                else: 
                    Tm[i] = Pm[i]*Tm[i+1]+Qm[i]
        #Contador
        cont = cont+1
        #Calculo do Erro
        erro = np.zeros(n_volumes)
        for i in range(0,n_volumes):
            erro[i] = ((Tm[i] - Tm_i[i])**2)**(1/2)
            Tm_i[i]=Tm[i]
            
        erro = np.max(erro)    
          
        
    fim = time.time() # Tempo final de convergência  
    
    #Cálculo da psoição dos centros dos volumes na geometria
    x=np.zeros(n_volumes)
    for i in range(0,n_volumes):
        x[i] = (L/n_volumes)/2 + ((L/n_volumes)/2)*i*2
        
    print(f'Tempo : {fim - inicio}')
    return Tm,cont,x
#TDMA(T_inf,T_b,k,D,L,h,alfa,delta_t,Tp_zero,n_volumes,criterio): 

Temperaturas, cont, x = TDMA(293,373,10,0.01,0.05,5,1*10**(-6),0.5,293,10,0.001)

print(f'Temperaturas: {Temperaturas} \n')
print(f'cont: {cont} \n')
print(f'x: {x} \n')
