%#######################################################################
%
%                       femur_icng Program
%
%           Calculates the intercondylar notch geometry given the Osirix
%           region-of-interest comma separated files.  The measurements
%           include the overall notch volume and the width at a four
%           different locations, denoted AA, MA , O and I.
%
%           Notes:  1. Each volume file must follow a specific naming convention
%                      consisting of the labels 'notchvol' & 'notchvol2'
%
%           12March11 - Daniel Sturnick
%
%#######################################################################
sthickT1a = 0.3125;
sthickT1s = 1.2;
smissing = 0;
iplt = true;
%
% Get knee files
%
d = dir('*NW_DW.csv');
d = char(d.name);
d = reorder(d);
%
% Let user select Files
%
ok =0;
while (ok==0)
     [filesel,ok] = listdlg('ListString',d,'Name', ...
     'ROI Data Files','PromptString', ...
     'Select a file for curvature plotting/measurements', ...
     'SelectionMode','multiple');
end
fnams = d(filesel,:);
fnams = reorder(fnams);
nfiles = size(fnams,1);
fg = figure;
%
% Create .csv data file
%
fid = fopen('femur_icng.csv','w');
fprintf(fid,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n', ...
    'Subj#','Knee#','Case/Cont','Inj/Uninj','Leg','Sex','2nd Inj','Postrecon', ...
    'AR','Study#','MRI Day','Age','Ht.(cm)','Wt. (kg)','Dleg','BW','NW_AA', ...
    'NW_I','NW_MA','NW_O','N_Vol','ACL_Vol','Ridge');
%
for k = 1:nfiles
%
% Get Notch Width File
%
    fnam1 = deblank(fnams(k,:));
%
% Prepare figure
%
figure(fg);
hold on;
%
% Get knee information
%
    idu = strfind(fnam1,'_');
    knum = fnam1(1:(idu(1)+1));
    while size(knum,2)<5
        knum = ['0' knum];
    end
    knum = deblank(knum);
    leg = strcmp(knum(5),'R');
%
% Get NW roi data
%   
    roi1 = rd_roi3(fnam1);
    roinams1 = lower(char(roi1.name));
    roinams1 = deblank(roinams1);
    idc = zeros(5,1);
    idc(1) = strmatchn('bw',roinams1,'exact');
    idc(2) = strmatchn('nw_aa',roinams1,'exact');
    idc(3) = strmatchn('nw_i',roinams1,'exact');
    idc(4) = strmatchn('nw_ma',roinams1,'exact');
    idc(5) = strmatchn('nw_o',roinams1,'exact');
    if any(idc==0)
        %keyboard;
        fprintf(fid,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n', ...
            ' ',knum,' ',' ',' ',' ',' ',' ',' ',' ',' ');
    else
    dat1 = roi1(idc(1)).data';
    xyzbw = dat1{:};
    p1 = plot3(xyzbw(:,1),xyzbw(:,2),xyzbw(:,3),'b-','LineWidth',1.5);
    dat2 = roi1(idc(2)).data';
    xyznwaa = dat2{:};
    p2 = plot3(xyznwaa(:,1),xyznwaa(:,2),xyznwaa(:,3),'c-','LineWidth',1.5);
    dat3 = roi1(idc(3)).data';
    xyznwi = dat3{:};
    p3 = plot3(xyznwi(:,1),xyznwi(:,2),xyznwi(:,3),'y-','LineWidth',1.5);
    dat4 = roi1(idc(4)).data';
    xyznwma = dat4{:};
    p4 = plot3(xyznwma(:,1),xyznwma(:,2),xyznwma(:,3),'m-','LineWidth',1.5);
    dat5 = roi1(idc(5)).data';
    xyznwo = dat5{:};
    p5 = plot3(xyznwo(:,1),xyznwo(:,2),xyznwo(:,3),'g-','LineWidth',1.5);
%
% Get notch widths
%
    bw = distance(xyzbw(1,:),xyzbw(2,:));
    nwaa = distance(xyznwaa(1,:),xyznwaa(2,:));
    nwi = distance(xyznwi(1,:),xyznwi(2,:));
    nwma = distance(xyznwma(1,:),xyznwma(2,:));
    nwo = distance(xyznwo(1,:),xyznwo(2,:));
%
% Get Notch Volume Data
%
    fnam2 = strrep(fnam1,'COROB_NW_DW.csv','AX_notchvol_DW.csv');
    if ~exist(fnam2,'file')
        fnam2 = [knum '_AX_notchvol2_SQ.csv'];
    end
%    
    if exist(fnam2,'file')
        roi2 = rd_roi3(fnam2);
        rnams = char(roi2.name);
        nroi = size(rnams,1);
        datv = roi2(1).data;
        datv = datv';
        if nroi>1
            for r = 2:nroi
                datn = roi2(r).data';
                datv = [datv; datn];
            end
        end
            
        nslice = size(datv,1);
        at1 = zeros(nslice,1);
        figure(fg);
        p6 = plt_datsl(datv,'k-');
        for s = 1:nslice
            xyz = datv{s};
            [pxyz,norm_vec,coeff,v,score,pexp,res,sse] = plane_fit(xyz(:,1),xyz(:,2),xyz(:,3));
            xyz = xyz*v;
            at1(s) = polyarea(xyz(:,1),xyz(:,2));
        end  
        nvol = sum(at1)*sthickT1a;
    else
        nvol = 9999;
    end
%
% Get ACL Volume Data
%
    fnam3 = ['acl data/' knum '_SAGACL_NF.csv'];
    if exist(fnam3,'file')    
        roi3 = rd_roi3(fnam3);
        datacl = roi3.data;
        datacl = datacl';
        p7 = plt_datsl(datacl,'r-',2);
        nslice2 = size(datacl,1);
        at2 = zeros(nslice2,1);
        for s = 1:nslice2
            xyz = datacl{s};
            [pxyz,norm_vec,coeff,v,score,pexp,res,sse] = plane_fit(xyz(:,1),xyz(:,2),xyz(:,3));
            xyz = xyz*v;
            at2(s) = polyarea(xyz(:,1),xyz(:,2));
        end  
        aclvol = sum(at2)*sthickT1s;
    else
        aclvol = 9999;
    end
%
% Get Notch Ridge Data
%
    fnamr = ['ridge data/' knum '_RIDGE_AS.csv'];
    if ~exist(fnamr,'file')
        fnamr = strrep(fnamr,'AS','DS');
    end
    roir = rd_roi3(fnamr);
    roinamsr = char(roir.name);
    idc = zeros(2,1);
    idc(1) = strmatch('WALL',roinamsr,'exact');
    idc(2) = strmatch('RIDGE',roinamsr,'exact');
    if any(idc==0)
        error('ERROR in icnotch_rdg: ROIs missing or labeled incorrectly.');
    end
%
% Get wall and ridge data in ordered cell arrays
%
    datW = roir(idc(1)).data';
    datR = roir(idc(2)).data';
    nslice = size(datW,1);
    if  ~(size(datR,1)==nslice)
        error('ERROR in icnotch_rdg: Corresponding number of wall and ridge ROIs do not exist');
    end
    mind = zeros(nslice,1);
    for s = 1:nslice
        xyzw = datW{s};
        for n=1:nslice
            xyzr= datR{n};
            if(smpln_tst(xyzr,xyzw,1e-8))
                mind(s)=n;
                break;
            end
        end
        if iplt
            plot3(xyzw(:,1),xyzw(:,2),xyzw(:,3),'bs','MarkerSize',4);
            plot3(xyzr(:,1),xyzr(:,2),xyzr(:,3),'k.','MarkerSize',8);
            axis equal;
        end
    end
    if any(mind==0)
           error('Slices are not lining up right!');
    end
    datR = datR(mind);
    [subj,cas,inj,sex] = subj_inf2(knum);
 %
 % Loop through slices
 %
    maxrdg = zeros(nslice,1);
    for s = 1:nslice
        % Fit wall to line
        xyzw = datW{s};
        [pxyz,rvec,r,pca,pexp,res] = line_fit(xyzw(:,1),xyzw(:,2),xyzw(:,3));
        xyzr = datR{s};
        npts = size(xyzr,1);
        xyzrht=zeros(npts,3); t = zeros(npts,1); rhts = zeros(npts,1);
        % Loop through ridge pts to get pt of max perpendicular distance
        for p = 1:npts
             [xyzrht(p,:),t(p)] = pt2line(pxyz,rvec',xyzr(p,:));
             rhts(p) = distance(xyzr(p,:),xyzrht(p,:));
        end
        [rht,ind] = max(rhts);
        rdgpk = xyzr(ind,:);
        rdgln = [xyzrht(ind,:);xyzr(ind,:)];
        wallln = repmat(pxyz,2,1)+[min([pca(:,1); t(ind)]); ...
                        max([pca(:,1); t(ind)])]*rvec';
        nvec = diff(rdgln);
        dvec = res*nvec';
        idn = find(dvec<0);
        if t(ind)>0
          idn2 = pca(idn,1)>0;
        else
          idn2 = pca(idn,1)<0;
        end
        idn = idn(idn2);
        dmin = sum(res(idn,:).^2,2);
        [dmin,idmin] = max(dmin);
        idmn = idn(idmin);
        oppln = [xyzw(idmn,:); pxyz+pca(idmn,1)*rvec'];
        if iplt
            figure(fg);
            plot3(rdgln(:,1),rdgln(:,2),rdgln(:,3),'b-','LineWidth',2);
            if leg
             text((rdgpk(1,1)-1),rdgpk(1,2),rdgpk(1,3),int2str(s),'FontSize',12, ...
                'FontWeight','bold','Color','k')
            else
             text((rdgpk(1,1)+1),rdgpk(1,2),rdgpk(1,3),int2str(s),'FontSize',12, ...
                'FontWeight','bold','Color','k')   
            end
            plot3(wallln(:,1),wallln(:,2),wallln(:,3),'k-','LineWidth',1);
            plot3(xyzr(ind,1),xyzr(ind,2),xyzr(ind,3),'mo', ...
                'LineWidth',2,'MarkerSize',4);
            plot3(oppln(:,1),oppln(:,2),oppln(:,3),'r.-','LineWidth',2);
            view(3);
            axis equal;
        end
        maxrdg(s) = rht+sqrt(dmin);
        %fprintf(fid1,'%g,"%s",%g,%g,%g,%i,%g\n',subj,knum,cas,inj,sex,s,maxrdg(s));
    end
    [ridge,ids] = max(maxrdg);
%
% Print data
%
    [subj,cas,inj,sex,sur,prev,ar,incomp,canum,conum,dem,mri,studynum,arf] = subj_inf2(knum);
    fprintf(fid,'%g,%s,%g,%g,%g,%g,%g,%g,%g,%s,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n',subj,knum,cas,inj, ...
        leg,sex,prev,sur,ar,studynum,mri,dem{1},dem{2},dem{3},dem{4},bw,nwaa,nwi,nwma,nwo,nvol,aclvol,ridge);
%
% Figure work
%
    figure(fg);
    view(3);
    title([knum ' Notch/ACL'],'FontSize',18,'FontWeight','bold','Interpreter','none');
    if ~((nvol==9999)||(aclvol==9999))
        le = legend([p1,p2,p3,p4,p5,p6,p7],'BW','NWAA','NWI','NWMA','NWO','Notch','ACL');
    elseif (~nvol==9999)&&(aclvol==9999)
        le = legend([p1,p2,p3,p4,p5,p6],'BW','NWAA','NWI','NWMA','NWO','Notch');
    else
        le = legend([p1,p2,p3,p4,p5],'BW','NWAA','NWI','NWMA','NWO');
    end
    pfile = [knum '_notchacl.ps'];
    print('-dpsc2',pfile);
    if k == 1
        print('-dpsc2','notchacl.ps');
    else
        print('-dpsc2','notchacl.ps','-append');
    end
    %pause;
    clf(fg);
    end
    
    clf(fg);
    %close(fg);
    
end