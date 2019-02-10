%{
Comparison between 2 groups of dancers

 Developers:  Andres Felipe Hernandez Marulanda
              Natali Olaya Mira
              Alejandro Miranda Hernandez

February 10
email:andres.hernandezma@upb.edu.co
%}
clc, clear all, close all
opc=menu('Choose selection ','Emg', 'Expert Emg', 'Beginer Emg','Expert vs Begginer')
switch opc
    case 1
      clc
        clear all
        fid=fopen('emg_healthy.dat')                 
        p1=fread(fid,50000,'short')
        p1=detrend(p1);
        plot(p1)
        title('EMG SIGNAL')
        Fs=1000; 
        Ts=1/Fs;
        nyq=Fs/2;
        L=1000;
        % %%%%%%%%%%%%%%%%%%%
        X=fft(p1); 
        % %%%%%%%%%%%%%%%%
        Y=abs(X/L);
        Y2=Y(1:L/2); 
        f=Fs*(1:L/2)/L; 
        % %%%%%%%%%%%%%%%%%%%
        subplot(3,1,1)
        plot(f,Y2) 
        title('Signal Fourier transform')
        % %%%%%%%%%%%%%%%%%%%%%%%
        wo=60/nyq;
        bw=wo/35;
        [b,a] = iirnotch(wo,bw,35);
        subplot(3,1,2)
        g=filter(b,a,p1);
        % %%%%%%%%%%%%%%%%%%%
         X1=fft(g);
        % %%%%%%%%%%%%%%%%%%%
         Yn=abs(X1/L); 
         Y2n=Yn(1:L/2); 
        % %%%%%%%%%%%%%%%%%
        plot(f,Y2n,'r')
        title('Notch Filter')
         [B1,A1] = butter(1,20/nyq,'high')
         [B2,A2] = butter(1,500/nyq,'low')
         i=filter(B1,A1,g);
         j=filter(B2,A2,i);
         subplot(3,1,3)
        % %%%%%%%%%%%%%%%%%%%
         X2=fft(j);
        % %%%%%%%%%%%%%%%%%%%
         Ynn=abs(X2/L); 
         Y2nn=Ynn(1:L/2); 
         f=Fs*(1:L/2)/L;
        % %%%%%%%%%%%%%%%%%
         plot(f,Y2nn,'r')
         title('Final Filtered Signal')
        % 
    case 2
        load('expert.mat');
        for d=1:size(datastart,1)
            canales(d,:) = [data(1,datastart(d):dataend(d))];
        end
        figure('Name','expert');
        subplot(4,1,1);
        plot(canales(3,:),'r');
        title('Adductor');
        xlabel('Time [s]');
        ylabel('Amplitude [cm]');
        legend('EMG');
        subplot(4,1,2);
        plot(canales(4,:),'m');
        title('Soleus');
        xlabel('Time [s]');
        ylabel('Amplitude [cm]');
        legend('EMG');
        subplot(4,1,3);
        plot(canales(5,:),'g');
        title('Internal Gastrocnemius');
        xlabel('Time [s]');
        ylabel('Amplitude [cm]');
        legend('EMG');
        subplot(4,1,4);
        plot(canales(6,:),'b');
        title('External Gastrocnemius');
        xlabel('Time [s]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        Fs = 1000;
        Ts = 1/Fs;
        L = 46400;
        t = (1:L)*Ts;
        f = Fs*(1:L/2)/L;

        abduc = fft(canales(3,:));
        Abduc = abs(abduc/L);
        Adductor = Abduc(1:L/2);

        figure('Name','expert');
        subplot(4,1,1);
        plot(f,Adductor,'r');
        title('Fourier Adductor');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        sol = fft(canales(4,:));
        Sol = abs(sol/L);
        Soleus = Sol(1:L/2);
        subplot(4,1,2);
        plot(f,Soleus,'m');
        title('Fourier Soleus');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG')

        gastri = fft(canales(5,:));
        Gastri = abs(gastri/L);
        Gastroi = Gastri(1:L/2);
        subplot(4,1,3);
        plot(f,Gastroi,'g');
        title('Fourier Internal Gastrocnemius');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        gastrd = fft(canales(6,:));
        Gastrd = abs(gastrd/L);
        Gastrod = Gastrd(1:L/2);
        subplot(4,1,4);
        plot(f,Gastrod,'b');
        title('Fourier External Gastrocnemius');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        for dd=1:size(canales,1)
            s=canales(dd,:);
            wn=[20,499]/(Fs/2);
            [b,a]=butter(5,wn);
            sf1=filtfilt(b,a,s);
            W0 = 60/(Fs/2);
            BW = W0/35;
            [b1,a1] = iirnotch(W0,BW);
            sf(dd,:) = filtfilt(b1,a1,sf1);
        end


        abduc1 = fft(sf(3,:));
        Abduc1 = abs(abduc1/L);
        Adductor1 = Abduc1(1:L/2);
        figure('Name','Beginner');
        plot(f,Adductor1,'r');
        title('Adductor Filtered out ');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        sol1 = fft(sf(4,:));
        Sol1 = abs(sol/L);
        Soleus1 = Sol1(1:L/2);
        figure('Name','Beginner');
        plot(f,Soleus1,'m');
        title('Soleus Filtered out ');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        gastri1 = fft(sf(5,:));
        Gastri1 = abs(gastri1/L);
        Gastroi1 = Gastri1(1:L/2);
        figure('Name','Beginner');
        plot(f,Gastroi1,'g');
        title('Internal Gastrocnemius Filtered out ');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        gastrd1 = fft(sf(6,:));
        Gastrd1 = abs(gastrd1/L);
        Gastrod1 = Gastrd1(1:L/2);
        figure('Name','Beginner');
        plot(f,Gastrod1,'b');
        title('External Gastrocnemius Filtered out ');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        figure('Name','expert');
        subplot(4,1,1);
        plot(f,Adductor1,'r');
        title('Adductor Filtered out ');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        subplot(4,1,2);
        plot(f,Soleus1,'m');
        title('Soleus Filtered out ');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        subplot(4,1,3);
        plot(f,Gastroi1,'g');
        title('Internal Gastrocnemius Filtered out ');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        subplot(4,1,4);
        plot(f,Gastrod1,'b');
        title('External Gastrocnemius Filtered out ');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        for d=1:size(sf,1)
            j=0;
            for i=1:10:length(sf)-20
                j=j+1;
                chanel(d,j)=rms(sf(d,i:i+20));
            end

        end

        figure('Name','expert');
        subplot(4,1,1);
        plot(chanel(3,:),'r');
        title('Adductor RMS');
        xlabel('Time [s]');
        ylabel('Amplitude [cm]');
        legend('EMG');
        subplot(4,1,2);
        plot(chanel(4,:),'m');
        title('Soleus RMS');
        xlabel('Time [s]');
        ylabel('Amplitude [cm]');
        legend('EMG');
        subplot(4,1,3);
        plot(chanel(5,:),'g');
        title('Internal Gastrocnemius RMS');
        xlabel('Time [s]');
        ylabel('Amplitude [cm]');
        legend('EMG');
        subplot(4,1,4);
        plot(chanel(6,:),'b');
        title('External Gastrocnemius RMS');
        xlabel('Time [s]');
        ylabel('Amplitude [cm]');
        legend('EMG');
    case 3
        load('Beginner.mat');
        for bd=1:size(datastart,1)
            bcanales(bd,:) = [data(1,datastart(bd):dataend(bd))];
        end
        figure('Name','Beginner');
        subplot(4,1,1);
        plot(bcanales(3,:),'r');
        title('Adductor');
        xlabel('Time [s]');
        ylabel('Amplitude [cm]');
        legend('EMG');
        subplot(4,1,2);
        plot(bcanales(4,:),'m');
        title('Soleus');
        xlabel('Time [s]');
        ylabel('Amplitude [cm]');
        legend('EMG');
        subplot(4,1,3);
        plot(bcanales(5,:),'g');
        title('Internal Gastrocnemius');
        xlabel('Time [s]');
        ylabel('Amplitude [cm]');
        legend('EMG');
        subplot(4,1,4);
        plot(bcanales(6,:),'b');
        title('External Gastrocnemius');
        xlabel('Time [s]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        bFs = 1000;
        bTs = 1/bFs;
        bL = 46400;
        bt = (1:bL)*bTs;
        bf = bFs*(1:bL/2)/bL;

        babduc = fft(bcanales(3,:));
        bAbduc = abs(babduc/bL);
        bAdductor = bAbduc(1:bL/2);

        figure('Name','Beginner');
        subplot(4,1,1);
        plot(bAdductor,'r');
        title('Fourier Adductor');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        bsol = fft(bcanales(4,:));
        bSol = abs(bsol/bL);
        bSoleus = bSol(1:bL/2);
        subplot(4,1,2);
        plot(bf,bSoleus,'m');
        title('Fourier Soleus');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG')

        bgastri = fft(bcanales(5,:));
        bGastri = abs(bgastri/bL);
        bGastroi = bGastri(1:bL/2);
        subplot(4,1,3);
        plot(bf,bGastroi,'g');
        title(' Fourier Internal Gastrocnemius');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        bgastrd = fft(bcanales(6,:));
        bGastrd = abs(bgastrd/bL);
        bGastrod = bGastrd(1:bL/2);
        subplot(4,1,4);
        plot(bGastrod,'b');
        title('Fourier External Gastrocnemius');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        for bdd=1:size(bcanales,1)
            bs=bcanales(bdd,:);
            bwn=[20,499]/(bFs/2);
            [bb,ba]=butter(5,bwn);
            bsf1=filtfilt(bb,ba,bs);
            bW0 = 60/(bFs/2);
            bBW = bW0/35;
            [bb1,ba1] = iirnotch(bW0,bBW);
            bsf(bdd,:) = filtfilt(bb1,ba1,bsf1);
        end

        babduc1 = fft(bsf(3,:));
        bAbduc1 = abs(babduc1/bL);
        bAdductor1 = Abduc1(1:bL/2);
        figure('Name','Beginner');
        plot(bf,bAdductor1,'r');
        title('Adductor Filtered out ');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        bsol1 = fft(bsf(4,:));
        bSol1 = abs(bsol/bL);
        bSoleus1 = bSol1(1:bL/2);
        figure('Name','Beginner');
        plot(bf,bSoleus1,'m');
        title('Soleus Filtered out ');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        bgastri1 = fft(bsf(5,:));
        bGastri1 = abs(bgastri1/bL);
        bGastroi1 = bGastri1(1:bL/2);
        figure('Name','Beginner');
        plot(bf,bGastroi1,'g');
        title('Internal Gastrocnemius Filtered out ');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        bgastrd1 = fft(bsf(6,:));
        bGastrd1 = abs(bgastrd1/bL);
        bGastrod1 = bGastrd1(1:bL/2);
        figure('Name','Beginner');
        plot(bf,bGastrod1,'b');
        title('External Gastrocnemius Filtered out ');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG');




        figure('Name','Beginner');
        subplot(4,1,1);
        plot(bf,bAdductor1,'r');
        title('Adductor Filtered out ');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        subplot(4,1,2);
        plot(bf,bSoleus1,'m');
        title('Soleus Filtered out ');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        bgastri1 = fft(bsf(5,:));
        bGastri1 = abs(bgastri1/bL);
        bGastroi1 = bGastri1(1:bL/2);
        subplot(4,1,3);
        plot(bf,bGastroi1,'g');
        title('Internal Gastrocnemius Filtered out ');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        subplot(4,1,4);
        plot(bf,bGastrod1,'b');
        title('External Gastrocnemius Filtered out ');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        for bd=1:size(bsf,1)
            j=0;
            for i=1:10:length(bsf)-20
                j=j+1;
                bchanel(bd,j)=rms(bsf(bd,i:i+20));
            end

        end

        figure('Name','Beginner');
        subplot(4,1,1);
        plot(bchanel(3,:),'r');
        title('Adductor RMS');
        xlabel('Time [s]');
        ylabel('Amplitude [cm]');
        legend('EMG');
        subplot(4,1,2);
        plot(bchanel(4,:),'m');
        title('Soleus RMS');
        xlabel('Time [s]');
        ylabel('Amplitude [cm]');
        legend('EMG');
        subplot(4,1,3);
        plot(bchanel(5,:),'g');
        title('Internal Gastrocnemius RMS');
        xlabel('Time [s]');
        ylabel('Amplitude [cm]');
        legend('EMG');
        subplot(4,1,4);
        plot(bchanel(6,:),'b');
        title('External Gastrocnemius RMS');
        xlabel('Time [s]');
        ylabel('Amplitude [cm]');
        legend('EMG');
    case 4
        load('expert.mat');
        for d=1:size(datastart,1)
            canales(d,:) = [data(1,datastart(d):dataend(d))];
        end
        figure('Name','expert');
        subplot(4,1,1);
        plot(canales(3,:),'r');
        title('Adductor');
        xlabel('Time [s]');
        ylabel('Amplitude [cm]');
        legend('EMG');
        subplot(4,1,2);
        plot(canales(4,:),'m');
        title('Soleus');
        xlabel('Time [s]');
        ylabel('Amplitude [cm]');
        legend('EMG');
        subplot(4,1,3);
        plot(canales(5,:),'g');
        title('Internal Gastrocnemius');
        xlabel('Time [s]');
        ylabel('Amplitude [cm]');
        legend('EMG');
        subplot(4,1,4);
        plot(canales(6,:),'b');
        title('External Gastrocnemius');
        xlabel('Time [s]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        Fs = 1000;
        Ts = 1/Fs;
        L = 46400;
        t = (1:L)*Ts;
        f = Fs*(1:L/2)/L;

        abduc = fft(canales(3,:));
        Abduc = abs(abduc/L);
        Adductor = Abduc(1:L/2);

        figure('Name','expert');
        subplot(4,1,1);
        plot(f,Adductor,'r');
        title('Fourier Adductor');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        sol = fft(canales(4,:));
        Sol = abs(sol/L);
        Soleus = Sol(1:L/2);
        subplot(4,1,2);
        plot(f,Soleus,'m');
        title('Fourier Soleus');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG')

        gastri = fft(canales(5,:));
        Gastri = abs(gastri/L);
        Gastroi = Gastri(1:L/2);
        subplot(4,1,3);
        plot(f,Gastroi,'g');
        title('Fourier Internal Gastrocnemius');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        gastrd = fft(canales(6,:));
        Gastrd = abs(gastrd/L);
        Gastrod = Gastrd(1:L/2);
        subplot(4,1,4);
        plot(f,Gastrod,'b');
        title('Fourier External Gastrocnemius');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        for dd=1:size(canales,1)
            s=canales(dd,:);
            wn=[20,499]/(Fs/2);
            [b,a]=butter(5,wn);
            sf1=filtfilt(b,a,s);
            W0 = 60/(Fs/2);
            BW = W0/35;
            [b1,a1] = iirnotch(W0,BW);
            sf(dd,:) = filtfilt(b1,a1,sf1);
        end


        abduc1 = fft(sf(3,:));
        Abduc1 = abs(abduc1/L);
        Adductor1 = Abduc1(1:L/2);
        figure('Name','Beginner');
        plot(f,Adductor1,'r');
        title('Adductor Filtered out ');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        sol1 = fft(sf(4,:));
        Sol1 = abs(sol/L);
        Soleus1 = Sol1(1:L/2);
        figure('Name','Beginner');
        plot(f,Soleus1,'m');
        title('Soleus Filtered out ');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        gastri1 = fft(sf(5,:));
        Gastri1 = abs(gastri1/L);
        Gastroi1 = Gastri1(1:L/2);
        figure('Name','Beginner');
        plot(f,Gastroi1,'g');
        title('Internal Gastrocnemius Filtered out ');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        gastrd1 = fft(sf(6,:));
        Gastrd1 = abs(gastrd1/L);
        Gastrod1 = Gastrd1(1:L/2);
        figure('Name','Beginner');
        plot(f,Gastrod1,'b');
        title('External Gastrocnemius Filtered out ');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        figure('Name','expert');
        subplot(4,1,1);
        plot(f,Adductor1,'r');
        title('Adductor Filtered out ');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        subplot(4,1,2);
        plot(f,Soleus1,'m');
        title('Soleus Filtered out ');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        subplot(4,1,3);
        plot(f,Gastroi1,'g');
        title('Internal Gastrocnemius Filtered out ');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        subplot(4,1,4);
        plot(f,Gastrod1,'b');
        title('External Gastrocnemius Filtered out ');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        for d=1:size(sf,1)
            j=0;
            for i=1:10:length(sf)-20
                j=j+1;
                chanel(d,j)=rms(sf(d,i:i+20));
            end

        end

        figure('Name','expert');
        subplot(4,1,1);
        plot(chanel(3,:),'r');
        title('Adductor RMS');
        xlabel('Time [s]');
        ylabel('Amplitude [cm]');
        legend('EMG');
        subplot(4,1,2);
        plot(chanel(4,:),'m');
        title('Soleus RMS');
        xlabel('Time [s]');
        ylabel('Amplitude [cm]');
        legend('EMG');
        subplot(4,1,3);
        plot(chanel(5,:),'g');
        title('Internal Gastrocnemius RMS');
        xlabel('Time [s]');
        ylabel('Amplitude [cm]');
        legend('EMG');
        subplot(4,1,4);
        plot(chanel(6,:),'b');
        title('External Gastrocnemius RMS');
        xlabel('Time [s]');
        ylabel('Amplitude [cm]');
        legend('EMG');
        %%
        %clear all
        load('Beginner.mat');
        for bd=1:size(datastart,1)
            bcanales(bd,:) = [data(1,datastart(bd):dataend(bd))];
        end
        figure('Name','Beginner');
        subplot(4,1,1);
        plot(bcanales(3,:),'r');
        title('Adductor');
        xlabel('Time [s]');
        ylabel('Amplitude [cm]');
        legend('EMG');
        subplot(4,1,2);
        plot(bcanales(4,:),'m');
        title('Soleus');
        xlabel('Time [s]');
        ylabel('Amplitude [cm]');
        legend('EMG');
        subplot(4,1,3);
        plot(bcanales(5,:),'g');
        title('Internal Gastrocnemius');
        xlabel('Time [s]');
        ylabel('Amplitude [cm]');
        legend('EMG');
        subplot(4,1,4);
        plot(bcanales(6,:),'b');
        title('External Gastrocnemius');
        xlabel('Time [s]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        bFs = 1000;
        bTs = 1/bFs;
        bL = 46400;
        bt = (1:bL)*bTs;
        bf = bFs*(1:bL/2)/bL;

        babduc = fft(bcanales(3,:));
        bAbduc = abs(babduc/bL);
        bAdductor = bAbduc(1:bL/2);

        figure('Name','Beginner');
        subplot(4,1,1);
        plot(bAdductor,'r');
        title('Fourier Adductor');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        bsol = fft(bcanales(4,:));
        bSol = abs(bsol/bL);
        bSoleus = bSol(1:bL/2);
        subplot(4,1,2);
        plot(bf,bSoleus,'m');
        title('Fourier Soleus');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG')

        bgastri = fft(bcanales(5,:));
        bGastri = abs(bgastri/bL);
        bGastroi = bGastri(1:bL/2);
        subplot(4,1,3);
        plot(bf,bGastroi,'g');
        title(' Fourier Internal Gastrocnemius');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        bgastrd = fft(bcanales(6,:));
        bGastrd = abs(bgastrd/bL);
        bGastrod = bGastrd(1:bL/2);
        subplot(4,1,4);
        plot(bGastrod,'b');
        title('Fourier External Gastrocnemius');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        for bdd=1:size(bcanales,1)
            bs=bcanales(bdd,:);
            bwn=[20,499]/(bFs/2);
            [bb,ba]=butter(5,bwn);
            bsf1=filtfilt(bb,ba,bs);
            bW0 = 60/(bFs/2);
            bBW = bW0/35;
            [bb1,ba1] = iirnotch(bW0,bBW);
            bsf(bdd,:) = filtfilt(bb1,ba1,bsf1);
        end

        babduc1 = fft(bsf(3,:));
        bAbduc1 = abs(babduc1/bL);
        bAdductor1 = Abduc1(1:bL/2);
        figure('Name','Beginner');
        plot(bf,bAdductor1,'r');
        title('Adductor Filtered out ');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        bsol1 = fft(bsf(4,:));
        bSol1 = abs(bsol/bL);
        bSoleus1 = bSol1(1:bL/2);
        figure('Name','Beginner');
        plot(bf,bSoleus1,'m');
        title('Soleus Filtered out ');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        bgastri1 = fft(bsf(5,:));
        bGastri1 = abs(bgastri1/bL);
        bGastroi1 = bGastri1(1:bL/2);
        figure('Name','Beginner');
        plot(bf,bGastroi1,'g');
        title('Internal Gastrocnemius Filtered out ');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        bgastrd1 = fft(bsf(6,:));
        bGastrd1 = abs(bgastrd1/bL);
        bGastrod1 = bGastrd1(1:bL/2);
        figure('Name','Beginner');
        plot(bf,bGastrod1,'b');
        title('External Gastrocnemius Filtered out ');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG');




        figure('Name','Beginner');
        subplot(4,1,1);
        plot(bf,bAdductor1,'r');
        title('Adductor Filtered out ');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        subplot(4,1,2);
        plot(bf,bSoleus1,'m');
        title('Soleus Filtered out ');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        bgastri1 = fft(bsf(5,:));
        bGastri1 = abs(bgastri1/bL);
        bGastroi1 = bGastri1(1:bL/2);
        subplot(4,1,3);
        plot(bf,bGastroi1,'g');
        title('Internal Gastrocnemius Filtered out ');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        subplot(4,1,4);
        plot(bf,bGastrod1,'b');
        title('External Gastrocnemius Filtered out ');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        for bd=1:size(bsf,1)
            j=0;
            for i=1:10:length(bsf)-20
                j=j+1;
                bchanel(bd,j)=rms(bsf(bd,i:i+20));
            end

        end

        figure('Name','Beginner');
        subplot(4,1,1);
        plot(bchanel(3,:),'r');
        title('Adductor RMS');
        xlabel('Time [s]');
        ylabel('Amplitude [cm]');
        legend('EMG');
        subplot(4,1,2);
        plot(bchanel(4,:),'m');
        title('Soleus RMS');
        xlabel('Time [s]');
        ylabel('Amplitude [cm]');
        legend('EMG');
        subplot(4,1,3);
        plot(bchanel(5,:),'g');
        title('Internal Gastrocnemius RMS');
        xlabel('Time [s]');
        ylabel('Amplitude [cm]');
        legend('EMG');
        subplot(4,1,4);
        plot(bchanel(6,:),'b');
        title('External Gastrocnemius RMS');
        xlabel('Time [s]');
        ylabel('Amplitude [cm]');
        legend('EMG');

        %%comparacion
        figure('Name','Comparison');
        [fi,co]=size(bchanel);
        subplot(4,1,1);
        plot(t(1:co),chanel(3,1:co),'b',bt(1:co),bchanel(3,:),'r');
        title('Adductor Expert-Beginner RMS');
        xlabel('Time [s]');
        ylabel('Amplitude [cm]');
        legend('Expert','Beginner');
        subplot(4,1,2);
        plot(t(1:co),chanel(4,1:co),'b',bt(1:co),bchanel(4,:),'m');
        title('Soleus Expert-Beginner RMS');
        xlabel('Time [s]');
        ylabel('Amplitude [cm]');
        legend('Expert','Beginner');
        subplot(4,1,3);
        plot(t(1:co),chanel(5,1:co),'b',bt(1:co),bchanel(5,:),'g');
        title('Internal Gastrocnemius Expert-Beginner RMS');
        xlabel('Time [s]');
        ylabel('Amplitude [cm]');
        legend('Expert','Beginner');
        subplot(4,1,4);
        plot(t(1:co),chanel(6,1:co),'b',bt(1:co),bchanel(6,:),'k');
        title('External Gastrocnemius Expert-Beginner RMS');
        xlabel('Time [s]');
        ylabel('Amplitude [cm]');
        legend('Expert','Beginner');

        for i=1:5
            comp(i,:)=chanel(i,1:co)-bchanel(i,:);
        end

        E3 = rms(chanel(3,:));
        E4 = rms(chanel(4,:));
        E5 = rms(chanel(5,:));
        E6 = rms(chanel(6,:));

        B3 = rms(bchanel(3,:));
        B4 = rms(bchanel(4,:));
        B5 = rms(bchanel(5,:));
        B6 = rms(bchanel(6,:));

        RMS = table([E3;E4;E5;E6],[B3;B4;B5;B6],'VariableNames',{'ExperT','Beginner'},'RowNames',{'RMS_Adductor','RMS_Soleus','RMS_Internal Gastrocnemius','RMS_External Gastrocnemius'})
end