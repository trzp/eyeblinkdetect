classdef EyeBlinkDetect
    properties
        srate;
        notch;
        notchB;
        notchA;
        bpB;
        bpA;
        meanNum;
        meanNumf;
        detectNum;
        threshold; %80
    end
    
    methods
        function Init(obj,srate,notch,threshold) %notch: fase/true
            obj.srate = srate;
            obj.meanNum = floor(obj.srate/50);    %目标降采样到50Hz
            obj.meanNumf = ones(1,obj.meanNum)./obj.meanNum;  
            obj.detectNum = ones(1,12);
            obj.threshold = threshold;
            
            obj.notch = notch
            if obj.notch
                %notch
                Fo = 50;
                Q = 15; %Q=35
                BW = (Fo/(srate/2))/Q;
                [obj.notchB,obj.notchA] = iircomb(srate/Fo,BW,'notch');  
            end
            
            fs = obj.srate/2;
            Wp = [1/fs 5/fs];
            Ws = [0.5/fs 10/fs];
            [N,Wn] = cheb1ord(Wp,Ws,3,40);
            [obj.bpB,obj.bpA] = cheby1(N,0.5,Wn);
        end
        
        function pos = detect(obj,signal,ana_length) %eoch决定截取signal尾部的一段信号进行分析，避免滤波产生的头部畸变
        if obj.notch
            signal = filtfilt(obj.notchB,obj.notchA,signal);
        end

        signal = filtfilt(obj.bpB,obj.bpA,signal);
        signal = signal(:,end-ana_length+1:end);
        signal = filter(obj.meanNumf,1,signal);
        signal = signal(:,floor(obj.meanNum/2):obj.meanNum,end);    %均值降采样
        [locs,pks] = findpeaks(signal);
        [locs1,pks1] = findpeaks(-1*signal);
        base = zeros(1,length(signal));
        base(locs) = pks;
        base(locs1) = pks1;
        base = filter(obj.detectNum,1,base);
        base = base(6:12:end);
        pos = find(baes > 80);
        end   
    end
end