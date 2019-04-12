%% 3D parametric lattices
close all
clear


unit_cells = 4; % go for a 3x3x3 system
amp = 0.8*pi/3;
l = 0.2;

t = -pi/2:0.5:(pi*(3/2));
a = pi;
b = pi/2;
c = pi;
f = 1; % period = 2*pi


figure
hold on
axis square

xlabel('U')
ylabel('V')
zlabel('W')

surface = [];
% B, B' in u, w
u = t;
w = -amp*sin(f*t);
v = 0*ones(size(t)) + amp*sin(f*t);
points = zeros(4*length(t), 3);

unit_offset = 2*pi;
for i = 1:length(t)
    points(4*i-3,:) = [u(i), v(i) + l, w(i) - l];
    points(4*i-2,:) = [u(i), v(i) + l, w(i) + l];
    points(4*i-1,:) = [u(i), v(i) - l, w(i) + l];
    points(4*i,:) = [u(i), v(i) - l, w(i) - l];
end
scatter3(points(:,1), points(:,2), points(:, 3), '.k')
surface = [surface; points];
v = a*ones(size(t)) - amp*sin(f*t);
for i = 1:length(t)
    points(4*i-3,:) = [u(i), v(i) + l, w(i) - l];
    points(4*i-2,:) = [u(i), v(i) + l, w(i) + l];
    points(4*i-1,:) = [u(i), v(i) - l, w(i) + l];
    points(4*i,:) = [u(i), v(i) - l, w(i) - l];
end
scatter3(points(:,1), points(:,2), points(:, 3), '.k')
surface = [surface; points];

u = t;
w = amp*sin(f*t)+pi;
v = 0*ones(size(t)) + amp*sin(f*t);
for i = 1:length(t)
    points(4*i-3,:) = [u(i), v(i) + l, w(i) - l];
    points(4*i-2,:) = [u(i), v(i) + l, w(i) + l];
    points(4*i-1,:) = [u(i), v(i) - l, w(i) + l];
    points(4*i,:) = [u(i), v(i) - l, w(i) - l];
end
scatter3(points(:,1), points(:,2), points(:, 3), '.k')
surface = [surface; points];
v = a*ones(size(t)) - amp*sin(f*t);
for i = 1:length(t)
    points(4*i-3,:) = [u(i), v(i) + l, w(i) - l];
    points(4*i-2,:) = [u(i), v(i) + l, w(i) + l];
    points(4*i-1,:) = [u(i), v(i) - l, w(i) + l];
    points(4*i,:) = [u(i), v(i) - l, w(i) - l];
end
scatter3(points(:,1), points(:,2), points(:, 3), '.k')
surface = [surface; points];

u = amp*sin(f*t);
w = t;
v = a*ones(size(t)) + amp*sin(f*t);
for i = 1:length(t)
    points(4*i-3,:) = [u(i) - l, v(i) + l, w(i)];
    points(4*i-2,:) = [u(i) + l, v(i) + l, w(i)];
    points(4*i-1,:) = [u(i) + l, v(i) - l, w(i)];
    points(4*i,:) = [u(i) - l, v(i) - l, w(i)];
end
scatter3(points(:,1), points(:,2), points(:, 3), '.k')
surface = [surface; points];


% C, C'
w = t;
v = -amp*sin(f*t);
u = 0*ones(size(t)) + amp*sin(f*t);
for i = 1:length(t)
    points(4*i-3,:) = [u(i) - l, v(i) + l, w(i)];
    points(4*i-2,:) = [u(i) + l, v(i) + l, w(i)];
    points(4*i-1,:) = [u(i) + l, v(i) - l, w(i)];
    points(4*i,:) = [u(i) - l, v(i) - l, w(i)];
end
scatter3(points(:,1), points(:,2), points(:, 3), '.k')
surface = [surface; points];
u = a*ones(size(t)) - amp*sin(f*t);
for i = 1:length(t)
    points(4*i-3,:) = [u(i) - l, v(i) + l, w(i)];
    points(4*i-2,:) = [u(i) + l, v(i) + l, w(i)];
    points(4*i-1,:) = [u(i) + l, v(i) - l, w(i)];
    points(4*i,:) = [u(i) - l, v(i) - l, w(i)];
end
scatter3(points(:,1), points(:,2), points(:, 3), '.k')
surface = [surface; points];


w = t;
v = amp*sin(f*t) + pi;
u = a*ones(size(t)) - amp*sin(f*t);
for i = 1:length(t)
    points(4*i-3,:) = [u(i) - l, v(i) + l, w(i)];
    points(4*i-2,:) = [u(i) + l, v(i) + l, w(i)];
    points(4*i-1,:) = [u(i) + l, v(i) - l, w(i)];
    points(4*i,:) = [u(i) - l, v(i) - l, w(i)];
end
scatter3(points(:,1), points(:,2), points(:, 3), '.k')
surface = [surface; points];

w = amp*sin(f*t);
v = t;
u = 0*ones(size(t)) - amp*sin(f*t);
for i = 1:length(t)
    points(4*i-3,:) = [u(i) - l, v(i), w(i) + l];
    points(4*i-2,:) = [u(i) + l, v(i), w(i) + l];
    points(4*i-1,:) = [u(i) + l, v(i), w(i) - l];
    points(4*i,:) = [u(i) - l, v(i), w(i) - l];
end
scatter3(points(:,1), points(:,2), points(:, 3), '.k')
surface = [surface; points];
u = a*ones(size(t)) + amp*sin(f*t);
for i = 1:length(t)
    points(4*i-3,:) = [u(i) - l, v(i), w(i) + l];
    points(4*i-2,:) = [u(i) + l, v(i), w(i) + l];
    points(4*i-1,:) = [u(i) + l, v(i), w(i) - l];
    points(4*i,:) = [u(i) - l, v(i), w(i) - l];
end
scatter3(points(:,1), points(:,2), points(:, 3), '.k')
surface = [surface; points];

w = -amp*sin(f*t)+pi;
v = t;
u = 0*ones(size(t)) - amp*sin(f*t);
for i = 1:length(t)
    points(4*i-3,:) = [u(i) - l, v(i), w(i) + l];
    points(4*i-2,:) = [u(i) + l, v(i), w(i) + l];
    points(4*i-1,:) = [u(i) + l, v(i), w(i) - l];
    points(4*i,:) = [u(i) - l, v(i), w(i) - l];
end
scatter3(points(:,1), points(:,2), points(:, 3), '.k')
surface = [surface; points];
u = a*ones(size(t)) + amp*sin(f*t);
for i = 1:length(t)
    points(4*i-3,:) = [u(i) - l, v(i), w(i) + l];
    points(4*i-2,:) = [u(i) + l, v(i), w(i) + l];
    points(4*i-1,:) = [u(i) + l, v(i), w(i) - l];
    points(4*i,:) = [u(i) - l, v(i), w(i) - l];
end
scatter3(points(:,1), points(:,2), points(:, 3), '.k')
surface = [surface; points];
%% Multiplying the Unit Cell:
for j = 1:unit_cells
    % U direction
    holder = surface(:,1) + unit_offset;
    surface = [surface; [holder, surface(:,2:3)]];
    % V direction
    holder = surface(:,2) + unit_offset;
    surface = [surface; [surface(:,1), holder, surface(:,3)]];
    % W direction
    holder = surface(:,3) + unit_offset;
    surface = [surface; [surface(:,1:2), holder]];
end

%% Scaling and offsetting
offset = min(surface);
surface = surface - offset;
scaling = max(surface);
surface = surface./scaling;

%% Constructing Surface
% Using alphaShape
figure
shp = alphaShape(surface(:,1), surface(:,2), surface(:,3));
shp.Alpha = 2*l/max(scaling);
plot(shp)
[triangles, vertices] = boundaryFacets(shp);
stlwrite('test0.stl', triangles, vertices)