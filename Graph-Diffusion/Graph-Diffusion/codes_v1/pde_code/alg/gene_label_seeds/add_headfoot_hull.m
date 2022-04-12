function headfoot_points=add_headfoot_hull(seeds_location,input_im)
% before convexhull add nodes around head and foot
% return points which include seeds_location and new points from head and
% foot

%% head
head = seeds_location(4,:);
neck = seeds_location(3,:);

disx = abs(neck(1)-head(1));
disy = abs(neck(2)-head(2));
distance = sqrt(disx^2+disy^2);
% distance = 2*distance/3;
% distance = 3*distance/4;

minx = head(1) - distance + 1;
miny = head(2) - distance + 1;
maxx = head(1) + distance + 1;
maxy = head(2) + distance + 1;

% imshow(input_im);
% hold on;
% plot(minx,miny,'ro');
% plot(maxx,maxy,'ro');
% hold off;

temp=[minx,miny];
headfoot_points=[seeds_location;temp ];
temp=[maxx,maxy];
headfoot_points=[headfoot_points;temp ];
temp=[minx,maxy];
headfoot_points=[headfoot_points;temp ];
temp=[maxx,miny];
headfoot_points=[headfoot_points;temp ];

%% foot
ankle_R = seeds_location(19,:);
foot_R = seeds_location(20,:);
ankle_L = seeds_location(15,:);
foot_L = seeds_location(16,:);


disRx=abs(ankle_R(1)-foot_R(1));
disRy=abs(ankle_R(2)-foot_R(2));
disR=3*sqrt(disRx^2+disRy^2)/4;

minRx=foot_R(1)-disR+1;
minRy=foot_R(2)-disR+1;
maxRx=foot_R(1)+disR+1;
maxRy=foot_R(2)+disR+1;

temp=[minRx,minRy];
headfoot_points=[headfoot_points;temp ];
temp=[maxRx,maxRy];
headfoot_points=[headfoot_points;temp ];
temp=[minRx,maxRy];
headfoot_points=[headfoot_points;temp ];
temp=[maxRx,minRy];
headfoot_points=[headfoot_points;temp ];

disLx=abs(ankle_L(1)-foot_L(1));
disLy=abs(ankle_L(2)-foot_L(2));
disL=3*sqrt(disLx^2+disLy^2)/4;

minLx=foot_L(1)-disL+1;
minLy=foot_L(2)-disL+1;
maxLx=foot_L(1)+disL+1;
maxLy=foot_L(2)+disL+1;

temp=[minLx,minLy];
headfoot_points=[headfoot_points;temp ];
temp=[maxLx,maxLy];
headfoot_points=[headfoot_points;temp ];
temp=[minLx,maxLy];
headfoot_points=[headfoot_points;temp ];
temp=[maxLx,minLy];
headfoot_points=[headfoot_points;temp ];

