%########################################################################
%
%               * icnotch_rdg Program *
%
%       Program to investigate the anteriomedial portion of the femoral
%       intercondylar notch.  Read's Osirix ROI .csv files.  Fits a line to
%       the notch wall, with an ROI label of 'WALL'.  Finds the greatest
%       perpendicular line to 
%
%       Notes
%
%       21Feb2011 * Daniel Sturnick
%
%########################################################################
%plot?
iplt = true;
%
% observers
%
obs = ['DS'];
%
% Get list of .csv files
%
d = dir(['*RIDGE2_' obs '.csv']);
nams = char(d.name);

fg = figure;

ok = 0;
while (ok==0)
     [filesel,ok] = listdlg('ListString',nams,'Name', ...
     'ROI Axial CSV Data Files','PromptString', ...
     'Select a file for plotting','SelectionMode','multiple');
end
fnams = nams(filesel,:);
nknees = size(fnams,1);
%
% Set up data file
%
if ispc
    fid = fopen('icnotch_rdg_max4.csv','wt');
    fid1 = fopen('icnotch_rdg4.csv','wt');
else
    fid = fopen('icnotch_rdg_max4.csv','w');
    fid1 = fopen('icnotch_rdg4.csv','w');
end
fprintf(fid,'%s,%s,%s,%s,%s,%s,%s\n','Subj#','KneeNum','Case/Cont','Injured?','Sex','Max Slice','Max Ridge');
fprintf(fid1,'%s,%s,%s,%s,%s,%s,%s\n','Subj#','KneeNum','Case/Cont','Injured','Sex','Slice#','Max Ridge');
%
% Loop through knees
%
for k = 1:nknees
%
% Get knee information
%
    knum = fnams(k,1:5);
    if strcmp(knum(5),'R')
     leg = 1;
    else
     leg = 0;
    end
    ktitle = [knum '_RIDGE_' obs];
%
% Set up figure
%
    if iplt
        figure(fg);
        hold on;
    end
%
% Get ROI's
%
    fnam = fnams(k,:);
    roi = rd_roi3(fnam);
    roinams = char(roi.name);
    idc = zeros(2,1);
    idc(1) = strmatch('WALL',roinams,'exact');
    idc(2) = strmatch('RIDGE',roinams,'exact');
    if any(idc==0)
        error('ERROR in icnotch_rdg: ROIs missing or labeled incorrectly.');
    end
%
% Get wall and ridge data in ordered cell arrays
%
    datW = roi(idc(1)).data';
    datR = roi(idc(2)).data';
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
        if iplt&&(s==4)
            %plot3(xyzw(:,1),xyzw(:,2),xyzw(:,3),'bs','MarkerSize',4);
            plot(xyzw(:,1),xyzw(:,2),'bs','MarkerSize',4);
            %plot3(xyzr(:,1),xyzr(:,2),xyzr(:,3),'k.','MarkerSize',8);
            plot3(xyzr(:,1),xyzr(:,2),xyzr(:,3),'k.','MarkerSize',8);
            axis equal;
        end
    end
    if any(mind==0)
           error('Slices are not lining up right!');
    end
    datR = datR(mind);
    [subj,cas,inj,sex] = subj_inf(knum);
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
        if iplt&&(s==4)
            plot3(rdgln(:,1),rdgln(:,2),rdgln(:,3),'b-','LineWidth',2);
            if leg
             text((rdgpk(1,1)-1),rdgpk(1,2),rdgpk(1,3),int2str(s),'FontSize',12, ...
                'FontWeight','bold','Color','k')
            else
             text((rdgpk(1,1)+1),rdgpk(1,2),rdgpk(1,3),int2str(s),'FontSize',12, ...
                'FontWeight','bold','Color','k')   
            end
            %plot3(wallln(:,1),wallln(:,2),wallln(:,3),'k-','LineWidth',1);
            plot(wallln(:,1),wallln(:,2),'k-','LineWidth',1);
            %plot3(xyzr(ind,1),xyzr(ind,2),xyzr(ind,3),'mo', ...
            %    'LineWidth',2,'MarkerSize',4);
            plot(xyzr(ind,1),xyzr(ind,2),'mo', ...
                'LineWidth',2,'MarkerSize',4);
            %plot3(oppln(:,1),oppln(:,2),oppln(:,3),'r.-','LineWidth',2);
            plot(oppln(:,1),oppln(:,2),'r.-','LineWidth',2);
            view(3);
            axis equal;
        end
        maxrdg(s) = rht+sqrt(dmin);
        fprintf(fid1,'%g,"%s",%g,%g,%g,%i,%g\n',subj,knum,cas,inj,sex,s,maxrdg(s));
        
    end
    [ridge,ids] = max(maxrdg);
    fprintf(fid,'%g,"%s",%g,%g,%g,%i,%g\n',subj,knum,cas,inj,sex,ids,ridge);
    if iplt
        title(ktitle,'FontSize', 26,'FontWeight','bold','Interpreter','none');
        xlabel('X');ylabel('Y');zlabel('Z');
    end
    pause;
    clf(fg);
end
%fclose(fid);