% Currell Berry Utkarsh Garg
% Assignment 6 part 2
% PBVM code -- much borrowed from getRangeAndBearing.m
% 4/15/2015

%%
image_width = 427;
focal_length_in_mm = 4;
ccd_width_in_mm = 4; % given 1/4" ccd width, http://www.dpreview.com/articles/8095816568/sensorsizes
focal_length = (image_width * focal_length_in_mm/ccd_width_in_mm)/1.65; % http://www.cs.cornell.edu/~snavely/bundler/focal.html
actual_distance = 0.07184; %m distance of centroid to bottom left

%%
im1 = iread('twoFeet45.jpg','grey','double')
im2 = iread('oneFoot00.jpg','grey','double')

c1 = icorner(im1, 'nfeat', 4)
c2 = icorner(im2, 'nfeat', 4)


%% compute image 1 stuff
%% now we extract coordinates of the top four features (corners) for each image
x1coords = [c1(1).u, c1(2).u, c1(3).u,  c1(4).u];
y1coords = [c1(1).v, c1(2).v, c1(3).v, c1(4).v];

all1coords = cat(2, x1coords', y1coords')
[sorted1, indices1] = sort(y1coords, 'descend');
centroid1 = [sum(x1coords)./4,
           	sum(y1coords)./4]';

sortedcoords1 = all1coords(indices1,:)

bottomleft1 = [];
bottomright1 = [];
if((sortedcoords1(1,1) < sortedcoords1(2,1)))
   bottomleft1 = sortedcoords1(1,:) 
   bottomright1 = sortedcoords1(2,:)
else
   bottomleft1 = sortedcoords1(2,:)
   bottomright1 = sortedcoords1(1,:)
end

dist1 = centroid1-bottomleft1;
image_distance1 = sqrt(sum(dist1.^2,2));
range1 = focal_length*actual_distance./image_distance1;

%% compute image 2 stuff
x2coords = [c2(1).u, c2(2).u, c2(3).u,  c2(4).u];
y2coords = [c2(1).v, c2(2).v, c2(3).v, c2(4).v];

all2coords = cat(2, x2coords', y2coords')
[sorted2, indices2] = sort(y2coords, 'descend');
centroid2 = [sum(x2coords)./4,
           	sum(y2coords)./4]';

sortedcoords2 = all2coords(indices2,:)

bottomleft2 = [];
bottomright2 = [];
if((sortedcoords2(1,1) < sortedcoords2(2,1)))
   bottomleft2 = sortedcoords2(1,:) 
   bottomright2 = sortedcoords2(2,:)
else
   bottomleft2 = sortedcoords2(2,:)
   bottomright2 = sortedcoords2(1,:)
end

dist2 = centroid2-bottomleft2;
image_distance2 = sqrt(sum(dist2.^2,2));
range2 = focal_length*actual_distance./image_distance2;


%% madness below
c1 = icorner(im1,'nfeat', 9)
c2 = icorner(im2,'nfeat', 9)
c1d = [c1.u ; c1.v]
c2d = [c2.u ; c2.v]
fmat = fmatrix(c1d,c2d)

K = [1211.2959 0          657.15924
	 0         1206.00512 403.17667
	 0         0           1       ];

emat = K'*fmat*K
[U,S,V] = svd(emat)

W = [0 -1 0
	 1  0 0
	 0  0 1];
Z = [ 0 1 0
	 -1 0 0
	 0 0 0];

R = U*W'*V' % note W' = inv(W)
ey = atan2(-R(3,1), sqrt(R(3,2).^2 + R(3,3).^2))


dst = sqrt(range1^2+range2^2-2*range1*range2*cos(ey))
translation = dst*3.28
disp('translation in feet:')
translation
disp('rotation in degrees:')
newBearing = acosd((range1.^2 + dst.^2 - range2.^2)/(2.*range1.*dst))

disp('motor values to take')
% note that 5.45 seconds goes 1 meter
% note that 1.2 seconds on 100 turns 90 degrees
timetoturn = newBearing*(1.2/90)
timetogostraight = dst*(1/5.45)
fprintf( '100, 0, %f\n', timetoturn)
fprintf( '0, 0, 1\n')
fprintf( '100, 100, %f\n', timetogostraight )
fprintf( '0, 0, 1\n')
 
%% The plotting only works for the particular case specified below
%% below, we generate the plot for the twoFoot45/onefoot00 case.
dst = translation
figure
hold on
initAngle = 45;
initDist = 2;
paperloc = [2.5 0];
targetloc = [2.5 1];
startloc = [2.5-sqrt(2),sqrt(2)];
scatter(startloc(1),startloc(2),'g');
scatter(paperloc(1),paperloc(2),'b');
scatter(targetloc(1),targetloc(2),'r');
correctLine = [startloc ; targetloc];
line(correctLine(:,1),correctLine(:,2));
side3 = sqrt(initDist^2+dst^2-2*initDist*dst*cosd(newBearing));
oangle = asind((sind(newBearing)/side3)*dst);
flipangle = 180-45;
inangle = 180 - 135-22;
xdisp = cosd(inangle)*dst;
ydisp = sind(inangle)*dst;
newloc = [startloc(1)+xdisp, startloc(2)-ydisp];
scatter(newloc(1),newloc(2),'k');
line2coords = [startloc ; newloc];
line(line2coords(:,1),line2coords(:,2));
xaxis(0,3);
yaxis(0,3);