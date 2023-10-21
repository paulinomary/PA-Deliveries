%% Chapter 9 - How are digital TV programs compressed to allow broadcasting?
%
% F. Marques, M. Menezes, J. Ruiz 2008
%
% This is a companion file to the book "Applied Signal Processing",
% by T.Dutoit and F. Marques, Springer 2008.
%
% It is supposed to be run cell-by-cell, using the cell mode of
% MATLAB 6 and later versions. Search for "what are cells" in the
% product help of MATLAB to know how to use them.
%
% This file uses the image toolbox of MATLAB.

%% Proof of concepts

% In the following Section we are going to illustrate the main concepts
% behind compression of image sequences. Image sequences present temporal
% as well as spatial redundancy. In order to correctly exploit the temporal
% redundancy between temporally neighbor images, the motion present in the
% scene has to be estimated. Once the motion is estimated, the information
% in the image used as reference is motion compensated to produce a first
% estimation of the image to be coded. Motion estimation is performed by
% dividing the image into non-overlapping square blocks and estimating the
% motion of each block independently.

iptsetpref('ImshowBorder', 'tight');format compact;


%% Macroblock processing

% Typically, for motion estimation, images are partitioned into blocks of
% 16x16 pixels , which will be referred to as macroblocks (see Fig. 9.7).
% The macroblock partition is hierarchical with respect to the block
% partition (see Chapter 8) and, therefore, every macroblock contains 4
% blocks. As in the block case, images are padded to allow an exact
% partition in terms of macroblocks. 

for i=1:3,
    table{i} = imread(sprintf('seqs/table/gray/table_%03d_g.bmp',i-1));
    figure,imshow(table{i});
    addgridtofigure(size(table{i}),[16 16]);
    title(['Table tennis #' int2str(i-1)]);
end;

% In Section 9.1, we have discussed that usual coding systems work in close
% loop; that is, the encoder uses decoded data (instead of the original
% one) to perform the prediction steps. This ensures that the transmitter
% and the receiver work in the same conditions. However, for the sake of
% simplicity, in the following Sections we are going to use original data
% in the prediction steps. Results in close loop will be presented in
% Section 9.2.9 when the complete coding system will be analyzed.
    

%% Block matching motion estimation

% Now, we are going to analyze how our basic motion compensation unit (that
% is, a generic 16x16 pixel macroblock) is processed. Among all the
% different techniques for estimating the motion associated to a macroblock
% with respect to a given reference image, the Block Matching algorithm has
% been shown to present very good coding properties.

% The Block Matching algorithm (BM) works independently for each macroblock
% of the current image and, for each macroblock, it looks for the best
% representation in the reference image, assuming that the whole macroblock
% has only undergone a translational movement (see Section 9.1.1).
% Therefore, the selected macroblock is placed at various positions in the
% reference image (those positions defined by the search area and the
% search strategy) and the pixel values overlapped by the displaced
% macroblock are compared (matched) with the original macroblock pixel
% values. The vector representing the displacement leading to the best
% match is the motion vector assigned to the selected macroblock. The final
% result of the BM algorithm is a set of displacements (motion vectors),
% each one associated to a different macroblock of the current image.

% As commented in the Section 9.1.1, several aspects have to be fixed in a
% concrete implementation of the BM; namely: (i) the metric that defines
% the best match; (ii) the search area in the reference image, which is
% related to the maximum allowed speed of objects in the scene; (iii) the
% quantization step to represent the coordinates of the motion vector,
% which is related to the accuracy of the motion representation; and (iv)
% the search strategy in the parameter space, which defines how many
% possible solutions (that is, different displacements) are analyzed in
% order to select the motion vector.

% The MATLAB function estimatemotion performs the BM of a given image with
% respect to a reference image. In it, the previous aspects are implemented
% as follows: (i) the selected metric is the absolute difference; (ii) the
% search area can be set by the fourth parameter (actually, this parameter
% does not directly set the search area but the maximum displacement
% allowed to the motion vector, so a value of [16 12] indicates that the
% motion vector can have values from (-16,-12) to (16,12)); (iii) the
% quantization step has been set to 1 pixel; and (iv) as search strategy
% both a full search and an nstep-search (see Fig. 9.4) within the search
% area are implemented.

t = 2;
[mvv mvh] = estimatemotion(table{t},table{t-1},[16 16],[15 15],'fullsearch');

% Fig. 9.8 presents the result of computing the BM between frames #1 and #0
% of the Table Tennis sequence. BM motion vectors are plotted in white over
% each of the macroblocks of frame #1 (motion vector field). Since an
% exhaustive search is performed, the estimated motion vectors are the
% optimum ones under the selected metric, search area and quantization
% step.

figure,plotmv(table{t},mvv,mvh,[16 16]);

% In order to illustrate how these vectors are computed, let us analyze the
% specific case of two different macroblocks. In Fig. 9.9, we show the
% information associated to the first of these macroblocks which is placed
% 5 blocks from the left and 2 from the top:

posr = 2; posc = 5;
f = table{t}*0.2;
bs = 16;
f(bs*(posr-1)+1:bs*(posr-1)+bs,bs*(posc-1)+1:bs*(posc-1)+bs) = ...
    table{t}(bs*(posr-1)+1:bs*(posr-1)+bs,bs*(posc-1)+1:bs*(posc-1)+bs);
figure,imshow(f);

% Fig. 9.10 shows a magnification of the selected macroblock:

b = table{t}(bs*(posr-1)+1:bs*(posr-1)+bs,bs*(posc-1)+1:bs*(posc-1)+bs); 
figure,imshow(kron(b,uint8(ones(8))));

% Fig. 9.11 shows the search area in the reference image (frame #0). The
% search area is centered at the same position of the macroblock under
% analysis and, in this case, the search area has been set to 46x46 pixels
% (parameter sa set to [15 15]). The macroblock from the current image
% (frame #1) is compared to all possible 16x16 subimages within the search
% area in the reference image (frame #0) (that is, 31x31 = 961 possible
% positions for each macroblock, see Section 9.1.1)

sa = [15 15];
f = table{t-1}*0.2;
f(bs*(posr-1)+1-sa(1):bs*(posr-1)+bs+sa(1),bs*(posc-1)+1-sa(2):bs*(posc-1)+bs+sa(2)) = ...
    table{t-1}(bs*(posr-1)+1-sa(1):bs*(posr-1)+bs+sa(1),bs*(posc-1)+1-sa(2):bs*(posc-1)+bs+sa(2));
figure,imshow(f);

% Fig. 9.12 shows a magnification of the search area:

br = table{t-1}(bs*(posr-1)+1-sa(1):bs*(posr-1)+bs+sa(1),bs*(posc-1)+1-sa(1):bs*(posc-1)+bs+sa(2)); 
figure,imshow(kron(br,uint8(ones(8))));


% Fig. 9.13 shows the error prediction surface; that is, the function
% containing the metric values of the comparison between the macroblock
% under analysis and all possible 16x16 subimages within the search area.
% Darker values of the error function (lower points in the 3D surface)
% correspond to lower error values between the macroblock and all possible
% subimages and, thus, better candidates to be the reference for the
% current macroblock. 

% function estimate_macroblock does the same as estimatemotion % but for
% only 1 macroblock 
[bcomp, bvv, bvh, errorf] = estimate_macroblock(b, table{t-1}, [posr posc], [16 16], [15 15],'fullsearch');bvv,bvh,

errorf_min = min(min(errorf));
errorf_max = max(max(errorf)) - errorf_min;
errorf_color = round(255.*(errorf - errorf_min)./errorf_max);

figure,mesh(-15:15,-15:15,errorf); %,errorf_color);
colormap('default'); %colormap(gray); 
colorbar; view(-20,22);
hold on; pcolor(-15:15,-15:15,errorf);
axis([-15 15 -15 15]);

% In this case, even if the error function presents several local minima,
% there is a global minimum (which correspond to the movement of the ball
% between frames #0 and #1) at position [0,8] (motion vector). That is, the
% best possible reference subimage (within the search area) for this
% macroblock is situated 8 pixels below and 0 pixels to the right of the
% original macroblock position. As the selected search strategy
% (full-search) performs an exhaustive search through the entire search
% area, it is guaranteed that the global minimum is always found. Other
% search strategies such as the nstep-search (See Section 9.1.1), which
% visits less positions in the search area, may be able to obtain similar
% (or even the same) results as the exhaustive search. For instance, for
% the case of the macroblock at position [2, 5], an nstep-search also finds
% the global minimum, as shown in Fig. 9.14:
   
[bcomp2, bvv2, bvh2, errorf2] = estimate_macroblock(b, table{t-1}, [posr posc], [16 16], [15 15],'nstep');bvv,bvh,
figure,surf(-15:15,-15:15,errorf2);
colormap('default'); 
colorbar;view(-20,22);
hold on; pcolor(-15:15,-15:15,errorf2);
axis([-15 15 -15 15]);


% In this example, both search strategies implemented in the MATLAB
% function estimatemotion (full-search and nstep-search) lead to the same
% result, a motion vector of [0,8]. Fig. 9.15 shows the selected reference
% subimage applying a motion vector of [0,8] within the search range.

f = table{t-1}*0.2;
f(bs*(posr-1)+1-sa(1):bs*(posr-1)+bs+sa(1),bs*(posc-1)+1-sa(2):bs*(posc-1)+bs+sa(2)) = ...
    0.45*table{t-1}(bs*(posr-1)+1-sa(1):bs*(posr-1)+bs+sa(1),bs*(posc-1)+1-sa(2):bs*(posc-1)+bs+sa(2));
f(bs*(posr-1)+1+bvv:bs*(posr-1)+bs+bvv,bs*(posc-1)+1+bvh:bs*(posc-1)+bs+bvh) = ...
    table{t-1}(bs*(posr-1)+1+bvv:bs*(posr-1)+bs+bvv,bs*(posc-1)+1+bvh:bs*(posc-1)+bs+bvh);
figure,imshow(f);

% Fig. 9.16 represents the magnification of the reference subimage and the
% corresponding compensated error (between the reference subimage at frame
% #0 and the macroblock at frame #1). 	

br = table{t-1}(bs*(posr-1)+1+bvv:bs*(posr-1)+bs+bvv,bs*(posc-1)+1+bvh:bs*(posc-1)+bs+bvh);
figure,imshow(kron(br,uint8(ones(8))));

e = double(b)-double(br);
figure,imshow_merror(kron(e,ones(8)),250);
sum(sum(e.*e)),


% Through this Chapter, in order to present error images, an offset of 128
% is added to all pixel values and, if necessary, they are clipped to the
% [0, 255] range. This way, zero error is represented by a 128 value

% A second example is analyzed to further illustrate the algorithm. Fig.
% 9.17 and Fig. 9.18 show the information associated to a macroblock from
% frame #1 placed at a different position. In this case the position of the
% macroblock under analysis is 6 macroblocks from the top and 6 from the
% left:

t = 2;
posr = 6; posc = 6;
f = table{t}*0.2;
bs = 16;
f(bs*(posr-1)+1:bs*(posr-1)+bs,bs*(posc-1)+1:bs*(posc-1)+bs) = ...
    table{t}(bs*(posr-1)+1:bs*(posr-1)+bs,bs*(posc-1)+1:bs*(posc-1)+bs);
figure,imshow(f);

% Fig. 9.18 shows a magnification of the macroblock under analysis:

b = table{t}(bs*(posr-1)+1:bs*(posr-1)+bs,bs*(posc-1)+1:bs*(posc-1)+bs); 
figure,imshow(kron(b,uint8(ones(8))));

% The search area at frame #0 corresponding to the macroblock under
% analysis can be seen in Fig. 9.19 and Fig. 9.20.

sa = [15 15];
f = table{t-1}*0.2;
f(bs*(posr-1)+1-sa(1):bs*(posr-1)+bs+sa(1),bs*(posc-1)+1-sa(2):bs*(posc-1)+bs+sa(2)) = ...
    table{t-1}(bs*(posr-1)+1-sa(1):bs*(posr-1)+bs+sa(1),bs*(posc-1)+1-sa(2):bs*(posc-1)+bs+sa(2));
figure,imshow(f);

br = table{t-1}(bs*(posr-1)+1-sa(1):bs*(posr-1)+bs+sa(1),bs*(posc-1)+1-sa(1):bs*(posc-1)+bs+sa(2)); 
figure,imshow(kron(br,uint8(ones(8))));

% As in the previous example, we show the error surface (see Fig. 9.21 and
% Fig. 9.22) corresponding to both full-search and nstep-search
% strategies[bcomp, bvv, bvh, errorf] = estimate_macroblock(b, table{t-1},
% [posr posc], [16 16], [15 15],'fullsearch');bvv,bvh,

figure,mesh(-15:15,-15:15,errorf);
colormap('default'); colorbar; view(-20,22);
hold on; pcolor(-15:15,-15:15,errorf);
axis([-15 15 -15 15]);

% In this example, the full-search obtains a motion vector of [0,-2] while
% the nstep-search obtains a motion vector of [-14,-1]. Note that, in this
% case, the nstep-search strategy does not lead to the optimum result. In
% particular, for the nstep-search, the initial step in the search is 8
% pixels. The nine solutions that are evaluated overlooked the global
% minimum and following steps of the algorithm get trapped in a distant
% local minimum. 

[bcomp2, bvv2, bvh2, errorf2] = estimate_macroblock(b, table{t-1}, [posr posc], [16 16], [15 15],'nstep');
bvv, %#ok<NOPTS>
bvh, %#ok<NOPTS>
figure,surf(-15:15,-15:15,errorf2);
colormap('gray'); colorbar;
view(-20,22);
hold on; pcolor(-15:15,-15:15,errorf2);
axis([-15 15 -15 15]);

% Fig. 9.23 shows the selected reference subimage leading to a motion
% vector of components [0,-2] (computed with the full-search approach).

f = table{t-1}*0.2;
f(bs*(posr-1)+1-sa(1):bs*(posr-1)+bs+sa(1),bs*(posc-1)+1-sa(2):bs*(posc-1)+bs+sa(2)) = ...
    0.45*table{t-1}(bs*(posr-1)+1-sa(1):bs*(posr-1)+bs+sa(1),bs*(posc-1)+1-sa(2):bs*(posc-1)+bs+sa(2));
f(bs*(posr-1)+1+bvv:bs*(posr-1)+bs+bvv,bs*(posc-1)+1+bvh:bs*(posc-1)+bs+bvh) = ...
    table{t-1}(bs*(posr-1)+1+bvv:bs*(posr-1)+bs+bvv,bs*(posc-1)+1+bvh:bs*(posc-1)+bs+bvh);
figure,imshow(f);

% Fig. 9.24 represents the magnification of the reference subimage and the
% corresponding compensated error. 

br = table{t-1}(bs*(posr-1)+1+bvv:bs*(posr-1)+bs+bvv,bs*(posc-1)+1+bvh:bs*(posc-1)+bs+bvh);
figure,imshow(kron(br,uint8(ones(8))));

e = double(b)-double(br);
figure,imshow_merror(kron(e,ones(8)),70);
sum(sum(e.*e)),


% In the case of the suboptimal strategy (nstep-search), Fig. 9.26 shows
% the selected reference subimage (leading to a motion vector of [-14,-1]
% and the corresponding compensated error. As the selected reference
% subimage is not the optimum one (within the search range), the final
% compensated error energy is greater in the case nstep-search.

f = table{t-1}*0.2;
f(bs*(posr-1)+1-sa(1):bs*(posr-1)+bs+sa(1),bs*(posc-1)+1-sa(2):bs*(posc-1)+bs+sa(2)) = ...
    0.45*table{t-1}(bs*(posr-1)+1-sa(1):bs*(posr-1)+bs+sa(1),bs*(posc-1)+1-sa(2):bs*(posc-1)+bs+sa(2));
f(bs*(posr-1)+1+bvv2:bs*(posr-1)+bs+bvv2,bs*(posc-1)+1+bvh2:bs*(posc-1)+bs+bvh2) = ...
    table{t-1}(bs*(posr-1)+1+bvv2:bs*(posr-1)+bs+bvv2,bs*(posc-1)+1+bvh2:bs*(posc-1)+bs+bvh2);
figure,imshow(f);

br = table{t-1}(bs*(posr-1)+1+bvv2:bs*(posr-1)+bs+bvv2,bs*(posc-1)+1+bvh2:bs*(posc-1)+bs+bvh2);
figure,imshow(kron(br,uint8(ones(8))));

e = double(b)-double(br);
figure,imshow_merror(kron(e,ones(8)),70);
sum(sum(e.*e)),


% The use of suboptimal search strategies can lead to the selection of
% reference subimages which are not the best representation of the
% macroblock under study. This translates into an increase of the
% prediction error and, therefore, a decrease in the coding performance.
% Note that, even in the case of using an exhaustive search, the resulting
% motion vectors may not represent the real motion in the scene. Vectors
% are computed for a fixed partition and obtained by minimization of a
% given metric in a specific neighbourhood. These constrains impose a
% scenario that may not be coincident with the real one. This problem is
% common in areas with very similar texture that may lead to different
% solutions with close values. Another situation that produces such kind of
% results is that of having two objects (or portions of objects) in the
% same macroblock undergoing different motions.



%% Motion compensation

% The estimation of the original frame #1 is obtained by motion
% compensating the reference image (in this example the frame #0) using the
% motion vectors from the BM algorithm, which created the so-called motion
% compensated image (see Fig. 9.27). In our implementation, the MATLAB
% function applymotion is responsible for compensating the reference image
% with the motion vectors computed using the previous function
% estimatemotion.

t=2;
[mvv mvh] = estimatemotion(table{t},table{t-1},[16 16],[15 15],'fullsearch');
tablec{t} = applymotion(table{t-1},mvv,mvh,[16 16]);
figure,imshow(tablec{t});
figure,imshow(table{t});

% Fig. 9.28.a shows the error image between the original and the motion
% compensated image (compensation error). 

error{t} = double(table{t})-double(tablec{t});
figure,imshow_merror(error{t},250);

% As it can be seen, the errors in the compensated image are located in the
% moving objects of the scene. However, as the computed motion vectors
% compensate somehow this motion (even if it is not translational) the
% final compensation error is relatively small.

Ee{t} = sum(sum(error{t}.*error{t}));
Ee{t},

% To have an idea of the improvement obtained by using motion compensation
% with BM, the direct difference between images #1 and #0 is presented in
% Fig. 9.28.b. This result can be interpreted as motion compensating image
% #0 with a motion vector field in which all motion vectors are [0,0]. It
% can be observed that the compensation error without motion compensation
% is larger than in the previous case, and the energy of the error would
% increase around 420%.

error_nocompen{t} = double(table{t})-double(table{t-1});
figure,imshow_merror(error_nocompen{t},250);
Eenc{t} = sum(sum(error_nocompen{t}.*error_nocompen{t})),
disp(['Energy increased with respect to compensation error = ' int2str(round((Eenc{t}-Ee{t})/Ee{t}*100)) '%']);


% Finally, Fig. 9.29 presents the motion compensated image (and the
% original image) which is obtained by estimating the motion vectors using
% the nstep-search sub-optimal strategy. In addition, Fig. 9.30 shows the
% corresponding compensation error. In that case, it can be seen than the
% motion compensated image is not able to estimate the original image as
% well as with the full search and, therefore, the energy of the error is
% greater than in the full search case. Also, some of the errors in the
% nstep-search strategy can be easily seen in the compensated image. For
% instance, the macroblock in position [6,6] (analyzed in Section 9.2.2)
% does not correctly reconstruct the racket or the player's hand due to the
% nstep-search strategy failing to obtain the optimum reference subimage in
% the search area.

t=2;
[mvv2 mvh2] = estimatemotion(table{t},table{t-1},[16 16],[15 15],'nstep');
tablec2{t} = applymotion(table{t-1},mvv2,mvh2,[16 16]);
figure,imshow(tablec2{t});
disp('motion compensated image (nstep)');

error2{t} = double(table{t})-double(tablec2{t});
figure,imshow_merror(error2{t},250);

Ee2{t} = sum(sum(error2{t}.*error2{t}));
Ee2{t}, %#ok<NOPTS>
disp(['Energy increased with respect to no motion = ' int2str(round((Eenc{t}-Ee2{t})/Ee2{t}*100)) '%']);

disp(['Energy increased with respect to exhaustive search = ' int2str(round((Ee{t}-Ee2{t})/Ee2{t}*100)) '%']);




%% Selection of search area

% Let us now analyze the effect of varying the size of the search area in
% the BM process. As we have previously commented, the search area is
% related to the maximum speed that we expect from the objects in the
% scene. In the case of the Table Tennis sequence the fastest object is the
% ball. We use frames #1 and #0 of the same sequence. Fig. 9.31 presents,
% on the top row, the motion compensated images using a search area with
% values 24x24, 32x32 and 46x46. On the second row the corresponding
% compensation errors are presented.

t = 2; 
n = 1;
for i=[4 8 15],
    sa = [2*i+16 2*i+16];
    [vv vh] = estimatemotion(table{t},table{t-1},[16 16],[i i],'fullsearch');
    mc = applymotion(table{t-1},vv,vh,[16 16]);
    %subplot(2,3,n);
    figure,imshow(mc);
    disp(['SA = ' mat2str(sa)]);
    e = double(table{t}) - double(mc);
    ee = sum(sum(e.^2));
    figure,imshow_merror(e,250); 
    disp(['Error energy = ' int2str(ee)]);
    n = n +1;
end;

% As it can be seen, a search area as small as 24x24 (corresponding to a
% maximum motion vector displacement of [±4, ±4]) does not allow correctly
% compensating the ball since its motion has displaced it outside the
% search area. This can be easily seen on the first column of Fig. 9.31
% where the ball does not appear in the motion compensated image (so the
% error in the area of the ball is very high). However, when using a search
% area of 32x32 or 46x46, the BM is capable of following the ball movement
% and the corresponding motion compensated images provide smaller errors in
% the area of the ball.

% Of course, the use of larger search areas implies a trade-off: the
% quality of the motion vectors is increased at the expenses of analyzing a
% larger amount of possible displacements and, therefore, increasing the
% computational load of the motion estimation step.


%% Selection of reference image

% Let us now analyze the impact of using a reference image different from
% the previous one. In order to do that, we will select frame #30 as the
% image to be estimated and successively use frame #28, #26, and #24 as
% reference images. The original images at those time instants can be seen
% in Fig. 9.32.

for i=25:2:31,
    table{i} = imread(sprintf('seqs/table/gray/table_%03d_g.bmp',i-1));
    figure,imshow(table{i}); 
    title(['original frame #' int2str(i)]);
end;

% Fig. 9.33 presents the compensated images and the associated compensation
% error images when using as reference image frame #24, #26 and #28,
% respectively. 

n = 1;
t = 31;
for p=29:-2:25,
    [mvv,mvh] = estimatemotion(table{t},table{p},[16 16],[15 15],'fullsearch');
    mc = applymotion(table{p},mvv,mvh,[16 16]);
    figure,imshow(mc);
    disp(['Reference #' int2str(p)]);
    e = double(table{t}) - double(mc);
    ee = sum(sum(e.^2));
    figure,imshow_merror(e,250); 
    disp(['Energy = ' int2str(ee)]);
end;

% As it can be seen, the further away the two images, the higher the energy
% of the compensation error. This is due to the fact that distant images
% are less correlated and, therefore, the prediction step cannot exploit in
% an efficient manner the temporal correlation. Moreover, there is an
% object, the player's arm at the bottom right corner, in frame #30 which
% is not present in frames #24 and #26. Therefore, when using these frames
% as references, this object does not have any referent in the reference
% image and the motion compensation leads to a compensation error with
% higher energy.


%% Backward motion estimation

% The problem of having an object in the current frame without any referent
% in the previous images happens every time an object appears in the scene.
% It is clear that, in order to obtain a referent for this object, we
% should look into the future frames; that is, use a non-causal prediction.
% At this stage, we are only going to analyze the possibility of using past
% as well as future frames as reference images. We postpone to Section
% 9.2.9 the illustration of how this non-causal prediction can actually be
% implemented in the coding context (see Sec-tion 9.1.2).

% Here, we analyze the possibility of estimating the current frame using as
% reference image the following frame in the sequence. If non-causal
% prediction can be used (that is, motion compensate future frames is
% possible), the information in past and future images can be used
% separately or combined in order to obtain a better motion compensated
% image. These concepts are illustrated in the next Figures.

% In the first row of Fig. 9.34, we present the current image (frame #30)
% with the motion vector fields with respect to the previous (Fig. 9.34.a)
% and the following (Fig. 9.34.b) frames. In the second row, we present the
% motion compensated images using these motion vector fields. Finally, in
% the last row, the associated compensation errors are shown.

for i=29:33,
    table{i} = imread(sprintf('seqs/table/gray/table_%03d_g.bmp',i-1));
end;

t=31;

% forward estimation
[fmvv,fmvh] = estimatemotion(table{t},table{t-1},[16 16],[15 15],'fullsearch');
fmc = applymotion(table{t-1},fmvv,fmvh,[16 16]);
fe = double(table{t}) - double(fmc);
fee = sum(sum(fe.^2));
 
% backward estimation
[bmvv,bmvh] = estimatemotion(table{t},table{t+1},[16 16],[15 15],'fullsearch');
bmc = applymotion(table{t+1},bmvv,bmvh,[16 16]);
be = double(table{t}) - double(bmc);
bee = sum(sum(be.^2));

% show images
figure,plotmv(table{t},fmvv,fmvh,[16 16]);
disp(['MV with Reference #' int2str(t-1)]);

figure,plotmv(table{t},bmvv,bmvh,[16 16]);
disp(['MV with Reference #' int2str(t+1)]);

figure,imshow(fmc);
disp('Forward compensated image');

figure,imshow(bmc); 
disp('Backward compensated image');

figure,imshow_merror(fe,250); 
disp(['forward error energy = ' int2str(fee)]);

figure,imshow_merror(be,250); 
disp(['backward error energy = ' int2str(bee)]);


% We try also with backward  but this time references #28 and #32 are used

% % forward estimation
% [fmvv2,fmvh2] = estimatemotion(table{t},table{t-2},[16 16],[15 15],'fullsearch');
% fmc2 = applymotion(table{t-2},fmvv2,fmvh2,[16 16]);
% fe2 = double(table{t}) - double(fmc2);
% fee2 = sum(sum(fe2.^2));
% figure,plotmv(table{t},fmvv2,fmvh2,[16 16]);
% disp(['MV with Reference #' int2str(t-2)]);
%  
% % backward estimation
% [bmvv2,bmvh2] = estimatemotion(table{t},table{t+2},[16 16],[15 15],'fullsearch');
% bmc2 = applymotion(table{t+2},bmvv2,bmvh2,[16 16]);
% be2 = double(table{t}) - double(bmc2);
% bee2 = sum(sum(be2.^2));
% figure,plotmv(table{t},bmvv2,bmvh2,[16 16]);
% disp(['MV with Reference #' int2str(t+2)]);
% 
% figure,imshow(fmc2);
% disp('Forward compensated image');
% figure,imshow(bmc2); 
% disp('Backward compensated image');
% 
% figure,imshow_merror(fe2,250); 
% disp(['forward error energy = ' int2str(fee2)]);
% figure,imshow_merror(be2,250); 
% disp(['backward error energy = ' int2str(bee2)]);


% From the example above it is clear than backward motion estimation can
% reduce the compensation error especially when new objects appear in the
% scene. In this example, this is the case of the player’s arm as well as a
% part of the background (note that the camera is zooming-out). However, in
% order to further reduce the compensation error, both motion compensated
% images (forward and backward) could be combined. Following this new idea,
% Fig. 9.35 shows the motion compensated image and the compensation error
% where, for each macroblock, the estimation is computed by linearly
% combining the two reference subimages with equal weights (w = 0.5).
    
bdmc = uint8(0.5 * (double(fmc) + double(bmc)));
bde = double(table{t}) - double(bdmc);
bdee = sum(sum(bde.^2));
figure,imshow(bdmc);
figure,imshow_merror(bde,250); 
disp(['forward error energy = ' int2str(bdee)]);


% In practice, the combination between forward and backward prediction is
% done in a macroblock level. For each macroblock, a decision (based on
% minimum error) has been made to select a forward, backward or combined
% reference. Fig. 9.36 shows the motion compensated image and the
% compensation error where each macroblock selects the best prediction:

[bimc,bie,bim] = biestimatemotion(table{t},table{t-1},table{t+1},[16 16],[15 15],'fullsearch');
figure,imshow(bimc); 
disp('Bi-directional compensated image');

biee = sum(sum(double(bie).^2));
figure,imshow_merror(bie,250);
disp(['bi-directional error energy = ' int2str(biee)]);

% For this example, combining both forward and backward predictions reduce
% the energy of the compensation error up to 64% and 32% with respect to
% only forward or backward estimation, respectively.

round((fee-biee)/fee*100),
round((bee-biee)/bee*100),

% In order to be able to see which decision has been taken for each
% macroblocks (among forward, backward or combined estimation) Fig. 9.37
% shows the decision for each macroblock encoded in gray values. Note that
% backward estimation (black macroblocks) has been selected in almost all
% the macroblocks of the perimeter of the image. This is due to the fact
% that, in these images, the camera is performing a zoom-out and,
% therefore, new information comes into the scene at each image. This way,
% a correct referent for the information in this zone in frame #30 has to
% be looked for in future frames (frame #31, in this case). Note, as well,
% that for the largest amount of macroblocks the best estimations are built
% up by combining forward or backward information (macroblocks in grey).
% Finally, there exist a few macroblocks which are better motion
% compensated using only in-formation from the past (macroblocks in grey).
	
kbim = kron(bim, ones(16));
figure,imshow((kbim-1)/4);
addgridtofigure(size(kbim),[16 16]);

%% Coding of the compensation error  

% In all previous examples, it can be observed that, although motion
% compensated images show good quality, there is still relevant information
% in the compensation error. The compensation error can be understood as a
% still image that presents some correlation among adjacent pixels. Even
% though that correlation is much lower than in original images, the
% compensation error information is coded using the same strategy used when
% coding still images (see Chapter 8). This way, the motion compensation
% error image is partitioned into blocks of 8x8 pixels and each block is
% separately transformed (using a DCT transformation), quantized and
% entropy coded.

% The same transform (DCT) and quantization table are used to code each
% compensation error block (as seen in the still image coding case in
% Chapter 8). Moreover, different k values can be used to obtain different
% qualities in the reconstructed compensation error. In our implementation,
% the MATLAB functions quantizedct8diff and iquantizedct8diff perform the
% (direct and inverse) quantization of the DCT coefficients for a
% compensation error block. As commented in Section 9.1.2, all values in
% the quantization table Q are set to 16.  

% To illustrate the transformation and quantization steps, Fig. 9.38
% presents the reconstruction of the compensation error when quantizing the
% DCT coefficients of its blocks using k = 1, k = 3 and k = 5. The
% compensation error corresponds to the block matching between frames #28
% and #30 (using an exhaustive search) of the Table Tennis sequence (see
% Fig. 9.33 right column).

for i=29:33,
    table{i} = imread(sprintf('seqs/table/gray/table_%03d_g.bmp',i-1));
end;

t = 31;
[mvv,mvh] = estimatemotion(table{t},table{t-2},[16 16],[15 15],'fullsearch');
mc = applymotion(table{t-2},mvv,mvh,[16 16]);
e = double(table{t}) - double(mc);
%Q = jpegsteps; Q(1,1) = Q(1,1)*2;
Q = 16*ones(size(jpegsteps));
edct = double(blkproc(e, [8 8], @dct2));
for k=[1 3 5 10],
    edctq = blkproc(edct, [8 8], @quantizedct8diff, k.*Q);
    edctqi = blkproc(edctq,[8 8],@iquantizedct8diff,k.*Q);
    er = blkproc(edctqi, [8 8], @idct2);
    figure,imshow_merror(er,250);
end;

% As it can be seen, increasing the values in the quantization tables
% reduces the information in the decoded compensation error and, so does
% the quality of the reconstructed image (obtained adding the compensated
% image and the compensation error). In this example, this effect can be
% observed in the better representation of the texture of the background in
% the first image.

% Fig. 9.39 further illustrates this concept presenting the reconstructed
% images. 

for k=[1 3 5 10],
    edctq = blkproc(edct, [8 8], @quantizedct8diff, k.*Q);
    edctqi = blkproc(edctq,[8 8],@iquantizedct8diff,k.*Q);
    er = blkproc(edctqi, [8 8], @idct2);
    tabler = double(mc) + er;
    tabler = uint8(limit(tabler,0,255));
    figure,imshow(tabler);
    disp(['k = ' int2str(k) '; psnr = ' num2str(mypsnr(table{t},tabler))]);
end;

%% Entropy coding

% As in the still image coding case, the quantized DCT coefficients
% asso-ciated to each block of the compensation error are entropy coded. In
% addition, the motion vectors associated with each macroblock of the image
% are encoded to form the final bitstream (see Section 9.1.2).

% For the DCT coefficients, the same strategy as in still image coding is
% used. DC coefficients are decorrelated and encoded separately while AC
% coefficients are preceded with a runlength encoder. However, in the case
% of coding the compensation error, a different set of predefined Huffman
% codes are used. The MATLAB file "codes/immotiondiffcodes.mat" includes
% these Huffman codes corresponding to different values of DC coefficients
% (in variable dcdcode) and AC coefficients (variables nzcode and vvcode).
% In our implementation, the MATLAB function encodeimage quantizes and
% encodes the corresponding DCT coefficients.

% Finally, motion vectors are also entropy coded. In our implementation,
% the MATLAB function encodemv is responsible for coding the motion
% vectors. This function imports the file "codes/mvcodes.mat" which
% includes the Huffman codes corresponding to the various motion vector
% values and returns the corresponding bitstream.

% The following code shows how motion vectors and compensation error are
% encoded into a single bitstream. Applying it to the motion vectors
% estimated between frames #28 and #30 of the Table Tennis se-quence and
% the corresponding compensation error results in:


t = 31;
[mvv,mvh] = estimatemotion(table{t},table{t-2},[16 16],...
    [15 15],'fullsearch');
mc = applymotion(table{t-2},mvv,mvh,[16 16]);
e = double(table{t}) - double(mc);
edct = double(blkproc(e, [8 8], @dct2));
mvcm = encodemv(mvv,mvh);
[diffcm,eQ] = encodeimage(e,Q,'motiondifference');
bitstream = strcat(mvcm, diffcm);
final_bits = size(bitstream,2),
 
% which yields the following figure of bits per pixel used in the
% representation of the decoded grey image:

bitspixel = final_bits / (size(table{t},1)*size(table{t},2)),


% In order to compare the final compressed frame with the original one, the
% produced bitstream must be decoded. The following code (using functions
% decodemv and decodeimage) shows how the bitstream can be decoded back
% into the reconstructed image. As entropy coding is a lossless process,
% the PSNR of the resulting reconstructed image is the same as in the
% previous section (see Fig. 9.39 with k = 1).

[mvv2, mvh2, bitstream] = decodemv(bitstream, size(table{t}),[16 16]);
[eQ2 bitstream] = decodeimage(bitstream, size(table{t}),Q,'motiondifference');
mc = applymotion(table{t-2},mvv,mvh,[16 16]);
tabler = uint8(limit(double(mc) + double(eQ2),0,255));
mypsnr(table{t},tabler),

            
%% Video coding

% In this section we are going to illustrate the entire video coding chain
% with different test sequences. In these experiments, we are going to work
% with color images. As in the still image coding case (see Chapter 8), the
% YCbCr color space is used and, to exploit the limitations of the human
% visual system, chrominance components are sub-sampled by a factor of two
% in the horizontal and vertical directions. The coding of chrominance
% components is done following the same steps as the luminance (grayscale)
% component. As it is common, the motion estimation is performed only in
% the luminance component and sub-sampled by two to be applied to the
% chrominance components.

% In our implementation, the MATLAB function mpegencode is responsible of
% performing the motion estimation, computing and coding the prediction
% error, closing the coding loop and selecting a frame type for the frame
% in the sequence. In order to force the use of I frames at the beginning
% of each GOP, the function is called separately for each GOP in the
% sequence. As this function does not implement any rate control, the
% quantization tables used in the coding of the compensation error have
% been multiplied by a factor of 2 to ensure that the resulting PSNR for P
% and B frames is similar to that of I frames.

% The following MATLAB code is used to encode the first 40 frames (5 GOPs
% of 8 frames) of the Table Tennis sequence. 5 different experiments are
% performed; using a GOP structure with no B frames (B=0), using 1 B frame
% between P frames (B=1), using 2 B frames between P frames (B=2), using 3
% B frames between P frames (B=3) and, finally, using 4 B frames between P
% frames (B=4). 

ngops = 5;
gopsize = 8;

for q=[1 2 3 4 5],
    
    for ng = 1:ngops,
    
        [cm,bits{1},type{1},psnr{1}] = mpegencode('results/table_color','seqs/table/color',1+(ng-1)*gopsize,ng*gopsize+1,  1,gopsize,q);
        [cm,bits{2},type{2},psnr{2}] = mpegencode('results/table_color','seqs/table/color',1+(ng-1)*gopsize,ng*gopsize+1,2,round(gopsize/2),q);
        [cm,bits{3},type{3},psnr{3}] = mpegencode('results/table_color','seqs/table/color',1+(ng-1)*gopsize,ng*gopsize+1,3,3,q);
        [cm,bits{4},type{4},psnr{4}] = mpegencode('results/table_color','seqs/table/color',1+(ng-1)*gopsize,ng*gopsize+1,4,2,q);
        [cm,bits{5},type{5},psnr{5}] = mpegencode('results/table_color','seqs/table/color',1+(ng-1)*gopsize,ng*gopsize+1,5,2,q);
        [cm,bits{6},type{6},psnr{6}] = mpegencode('results/table_color','seqs/table/color',1+(ng-1)*gopsize,ng*gopsize+1,6,2,q);

        for i=1:6,
            table{q}.bits{i}(1+(ng-1)*gopsize:ng*gopsize) = bits{i}(1+(ng-1)*gopsize:ng*gopsize);
            table{q}.type{i}(1+(ng-1)*gopsize:ng*gopsize) = type{i}(1+(ng-1)*gopsize:ng*gopsize);
            table{q}.psnr{i}.y(1+(ng-1)*gopsize:ng*gopsize) = psnr{i}.y(1+(ng-1)*gopsize:ng*gopsize);
            table{q}.psnr{i}.cb(1+(ng-1)*gopsize:ng*gopsize) = psnr{i}.cb(1+(ng-1)*gopsize:ng*gopsize);
            table{q}.psnr{i}.cr(1+(ng-1)*gopsize:ng*gopsize) = psnr{i}.cr(1+(ng-1)*gopsize:ng*gopsize);
        end;

        if (ng==ngops),
            for i=1:6,
                table{q}.bits{i}(ng*gopsize+1) = bits{i}(ng*gopsize+1);
                table{q}.type{i}(ng*gopsize+1) = type{i}(ng*gopsize+1);
                table{q}.psnr{i}.y(ng*gopsize+1) = psnr{i}.y(ng*gopsize+1);
                table{q}.psnr{i}.cb(ng*gopsize+1) = psnr{i}.cb(ng*gopsize+1);
                table{q}.psnr{i}.cr(ng*gopsize+1) = psnr{i}.cr(ng*gopsize+1);
            end;
        end;
    end;
    
end;


% Fig. 9.40 and Fig. 9.41 show the results of the 5 experiments for the
% sequence Table Tennis. In Fig. 9.41, the number of bits needed to encode
% each frame (for all 5 experiments) is represented and the corresponding
% PSNR (for the luminance component) is shown in Fig. 9.41. Furthermore,
% the average number of bits needed to represent the 40 frames and the
% average PSNR for all 40 frames are shown in Table 9.1.

% Let us discuss the coding efficiency concept by analyzing Fig. 9.40 and
% Fig. 9.41. Coding efficiency is measured in terms of achieving the same
% or better quality with the same or lower bitrate. First, note that, as
% already commented in Section 9.1.2, I-Pictures require many more bits
% than P-Pictures and B-Pictures, while achieving similar PSNR (the largest
% difference being 2 dBs). A similar conclusion can be driven with respect
% to P-Pictures and B-Pictures: due to the better compensation obtained in
% B-Pictures, slightly better qualities can be obtained with a lower number
% of bits. Note, as well, that when the number of B frames in the GOP
% increases, so does the number of bits necessary to encode the following
% P-Frame (for roughly the same PSNR). As we have discussed in Section
% 9.2.5, a larger number of B-Pictures translates into a P-Image which is
% more decorrelated with respect to its reference and, therefore, more bits
% are necessary to code the compensation error. Nevertheless, for the Table
% Tennis sequence, it can be seen that using B frames in the GOP structure
% really improves the coding efficiency of the video codec. 

% Finally, another aspect that can be commented in Fig. 9.41 is the
% degradation that the prediction chain introduces. This can be very well
% observed in the curve associated to B = 0 where the PSNR decays at each
% new predicted frame, even though the number of bits is progressively
% increased. 


% Fig. 9.42 represents the average PSNR and kbps needed to encode the 40
% frames of the Table Tennis sequence for each GOP structure and at
% different qualities. Qualities are selected by using a k factor which
% multiplies all values in the quantization tables (both for intra and
% inter prediction). 

% Fig. 9.42 shows that, for this sequence, using B frames generally
% improves the coding efficiency at high/normal qualities. However, for
% lower qualities, it is more efficient to encode using a lower number of B
% frames (ideally 0). This is due to the fact that, at lower qualities,
% reference frames in the closed loop are of very low quality and therefore
% the motion compensated images are not very similar to the frame being
% encoded. In this case, it is more efficient to use reference frames as
% close as possible to the frame being coded.

q = 1;
fcolor1 = {'bx-','ko-','rd-','g*-','c+-','ms-'};
fcolor2 = {'b','k','r','g','c','m'};
figure,grid,hold on;
for i=1:1:5,
    plot(table{q}.bits{i}/1000,fcolor1{i},'linewidth',1); 
    legend_str{i} = ['B = ' num2str(i-1)];  
end;
legend(legend_str);
for i=1:1:5,
%    plot(mean(table{q}.bits{i})*ones(size(table{q}.bits{i}))/1000,fcolor2{i},'linewidth',3);
    disp(['mean (' legend_str{i} ') = ' num2str(mean(table{q}.bits{i})/1000) ' [kbits]']);
end;
xlabel('Frame #'); ylabel('kbits');

% plot psnr2
figure,grid,hold on;
for i=1:1:5,
    plot(table{q}.psnr{i}.y,fcolor1{i},'linewidth',1); 
end;
legend(legend_str,'Location','NorthWest');
for i=1:1:5,
%    plot(meanpsnr(table{q}.psnr{i}.y)*ones(size(table{q}.psnr{i}.y)),fcolor2{i},'linewidth',3);
    disp(['meanpsnr (' legend_str{i} ') = ' num2str(meanpsnr(table{q}.psnr{i}.y)) ' [dB]']);
end;
xlabel('Frame #'); ylabel('PSNR_Y [dB]');

figure,grid,hold on;
for i=1:5,
    p = 1;
    for q=[1 2 3 4 5],
        val_bits(p) = sum(table{q}.bits{i}/1000)/length(table{q}.bits{i})*30;
        val_psnr(p) = meanpsnr(table{q}.psnr{i}.y);
        p = p+1;
    end;
    plot(val_bits,val_psnr,fcolor1{i},'MarkerSize',8,'LineWidth',2.2);
end;
legend(legend_str,'Location','SouthEast');
xlabel('kbps'); ylabel('PSNR_Y [dB]');


% Fig. 9.43, Fig. 9.44 and Fig. 9.45 show the same experiments con-ducted
% above but, this time, for 40 frames of the Foreman sequence. In this
% case, a specific part of the sequence is selected where there is a high
% amount of motion present. Similar conclusions can be driven as in the
% previous example. Nevertheless, as the motion estimation/compensation is
% not able to reduce the compensation error as much as in the previous test
% sequence, the more efficient GOP structure (at high/medium qualities)
% consists in using only 1 B frame between P frames (B = 1). Fig. 9.45
% shows (as in the previous test sequence) that, at lower qualities, the
% use of B frames reduces its coding efficiency. As previously, the average
% number of bits needed to represent the 40 frames and the average PSNR for
% all 40 frames are shown in Table 9.2.


ngops = 5;
gopsize = 8;

for q = [1 2 3 4 5],

    for ng = 1:ngops,
    
        [cm,bits{1},type{1},psnr{1}] = mpegencode('results/foreman_color','seqs/foreman/color_qcif',175+1+(ng-1)*gopsize,175+ng*gopsize+1,  1,gopsize,q);
        [cm,bits{2},type{2},psnr{2}] = mpegencode('results/foreman_color','seqs/foreman/color_qcif',175+1+(ng-1)*gopsize,175+ng*gopsize+1,2,round(gopsize/2),q);
        [cm,bits{3},type{3},psnr{3}] = mpegencode('results/foreman_color','seqs/foreman/color_qcif',175+1+(ng-1)*gopsize,175+ng*gopsize+1,3,3,q);
        [cm,bits{4},type{4},psnr{4}] = mpegencode('results/foreman_color','seqs/foreman/color_qcif',175+1+(ng-1)*gopsize,175+ng*gopsize+1,4,2,q);
        [cm,bits{5},type{5},psnr{5}] = mpegencode('results/foreman_color','seqs/foreman/color_qcif',175+1+(ng-1)*gopsize,175+ng*gopsize+1,5,2,q);
        [cm,bits{6},type{6},psnr{6}] = mpegencode('results/foreman_color','seqs/foreman/color_qcif',175+1+(ng-1)*gopsize,175+ng*gopsize+1,6,2,q);

        for i=1:6,
            foreman{q}.bits{i}(1+(ng-1)*gopsize:ng*gopsize) = bits{i}(175+1+(ng-1)*gopsize:175+ng*gopsize);
            foreman{q}.type{i}(1+(ng-1)*gopsize:ng*gopsize) = type{i}(175+1+(ng-1)*gopsize:175+ng*gopsize);
            foreman{q}.psnr{i}.y(1+(ng-1)*gopsize:ng*gopsize) = psnr{i}.y(175+1+(ng-1)*gopsize:175+ng*gopsize);
            foreman{q}.psnr{i}.cb(1+(ng-1)*gopsize:ng*gopsize) = psnr{i}.cb(175+1+(ng-1)*gopsize:175+ng*gopsize);
            foreman{q}.psnr{i}.cr(1+(ng-1)*gopsize:ng*gopsize) = psnr{i}.cr(175+1+(ng-1)*gopsize:175+ng*gopsize);
        end;

        if (ng==ngops),
            for i=1:6,
                foreman{q}.bits{i}(ng*gopsize+1) = bits{i}(175+ng*gopsize+1);
                foreman{q}.type{i}(ng*gopsize+1) = type{i}(175+ng*gopsize+1);
                foreman{q}.psnr{i}.y(ng*gopsize+1) = psnr{i}.y(175+ng*gopsize+1);
                foreman{q}.psnr{i}.cb(ng*gopsize+1) = psnr{i}.cb(175+ng*gopsize+1);
                foreman{q}.psnr{i}.cr(ng*gopsize+1) = psnr{i}.cr(175+ng*gopsize+1);
            end;
        end;

    end;
    
end;


%% plot bits (foreman)
q = 1;
fcolor1 = {'bx-','ko-','rd-','g*-','c+-','ms-'};
fcolor2 = {'b','k','r','g','c','m'};
figure,grid,hold on;
for i=1:1:5,
    plot(foreman{q}.bits{i}/1000,fcolor1{i},'linewidth',1); 
    legend_str{i} = ['B = ' num2str(i-1)];  
end;
legend(legend_str);
for i=1:1:5,
%    plot(mean(foreman{q}.bits{i})*ones(size(foreman{q}.bits{i}))/1000,fcolor2{i},'linewidth',3);
    disp(['mean (' legend_str{i} ') = ' num2str(mean(foreman{q}.bits{i})/1000) ' [kbits]']);
end;
xlabel('Frame #'); ylabel('kbits');
    
% plot psnr2
figure,grid,hold on;
for i=1:1:5,
    plot(foreman{q}.psnr{i}.y,fcolor1{i},'linewidth',1); 
end;
legend(legend_str,'Location','NorthWest');
for i=1:1:5,
%    plot(meanpsnr(foreman{q}.psnr{i}.y)*ones(size(foreman{q}.psnr{i}.y)),fcolor2{i},'linewidth',3);
    disp(['meanpsnr (' legend_str{i} ') = ' num2str(meanpsnr(foreman{q}.psnr{i}.y)) ' [dB]']);
end;
xlabel('Frame #'); ylabel('PSNR_Y [dB]');

figure,grid,hold on;
for i=1:5,
    p = 1;
    for q=[1 2 3 4 5],
        val_bits(p) = sum(foreman{q}.bits{i}/1000)/length(foreman{q}.bits{i})*30;
        val_psnr(p) = meanpsnr(foreman{q}.psnr{i}.y);
        p = p+1;
    end;
    plot(val_bits,val_psnr,fcolor1{i},'MarkerSize',8,'LineWidth',2.2);
end;
legend(legend_str,'Location','SouthEast');
xlabel('kbps'); ylabel('PSNR_Y [dB]');

ngops = 5;
gopsize = 8;

for ng = 1:ngops,
    
    [cm,bits{1},type{1},psnr{1}] = mpegencode('results/akiyo_color','seqs/akiyo/color_qcif',1+(ng-1)*gopsize,ng*gopsize+1,  1,gopsize);
    [cm,bits{2},type{2},psnr{2}] = mpegencode('results/akiyo_color','seqs/akiyo/color_qcif',1+(ng-1)*gopsize,ng*gopsize+1,2,round(gopsize/2));
    [cm,bits{3},type{3},psnr{3}] = mpegencode('results/akiyo_color','seqs/akiyo/color_qcif',1+(ng-1)*gopsize,ng*gopsize+1,3,3);
    [cm,bits{4},type{4},psnr{4}] = mpegencode('results/akiyo_color','seqs/akiyo/color_qcif',1+(ng-1)*gopsize,ng*gopsize+1,4,2);
    [cm,bits{5},type{5},psnr{5}] = mpegencode('results/akiyo_color','seqs/akiyo/color_qcif',1+(ng-1)*gopsize,ng*gopsize+1,5,2);
    [cm,bits{6},type{6},psnr{6}] = mpegencode('results/akiyo_color','seqs/akiyo/color_qcif',1+(ng-1)*gopsize,ng*gopsize+1,6,2);

    for i=1:6,
        akiyo.bits{i}(1+(ng-1)*gopsize:ng*gopsize) = bits{i}(1+(ng-1)*gopsize:ng*gopsize);
        akiyo.type{i}(1+(ng-1)*gopsize:ng*gopsize) = type{i}(1+(ng-1)*gopsize:ng*gopsize);
        akiyo.psnr{i}.y(1+(ng-1)*gopsize:ng*gopsize) = psnr{i}.y(1+(ng-1)*gopsize:ng*gopsize);
        akiyo.psnr{i}.cb(1+(ng-1)*gopsize:ng*gopsize) = psnr{i}.cb(1+(ng-1)*gopsize:ng*gopsize);
        akiyo.psnr{i}.cr(1+(ng-1)*gopsize:ng*gopsize) = psnr{i}.cr(1+(ng-1)*gopsize:ng*gopsize);
    end;
    
    if (ng==ngops),
        for i=1:6,
            akiyo.bits{i}(ng*gopsize+1) = bits{i}(ng*gopsize+1);
            akiyo.type{i}(ng*gopsize+1) = type{i}(ng*gopsize+1);
            akiyo.psnr{i}.y(ng*gopsize+1) = psnr{i}.y(ng*gopsize+1);
            akiyo.psnr{i}.cb(ng*gopsize+1) = psnr{i}.cb(ng*gopsize+1);
            akiyo.psnr{i}.cr(ng*gopsize+1) = psnr{i}.cr(ng*gopsize+1);
        end;
    end;
    
end;


%% plot bits (akiyo)
fcolor1 = {'bx-','ko-','rd-','g*-','c+-','ms-'};
fcolor2 = {'b','k','r','g','c','m'};
figure,grid,hold on;
for i=1:1:5,
    plot(akiyo.bits{i}/1000,fcolor1{i},'linewidth',1); 
    legend_str{i} = ['B = ' num2str(i-1)];  
end;
legend(legend_str);
for i=1:1:5,
    plot(mean(akiyo.bits{i})*ones(size(akiyo.bits{i}))/1000,fcolor2{i},'linewidth',3);
end;

xlabel('Frame #'); ylabel('kbits');
    

% plot psnr2
figure,grid,hold on;
for i=1:1:5,
    plot(akiyo.psnr{i}.y,fcolor1{i},'linewidth',1); 
end;
legend(legend_str,'Location','NorthWest');
for i=1:1:5,
    plot(meanpsnr(akiyo.psnr{i}.y)*ones(size(akiyo.psnr{i}.y)),fcolor2{i},'linewidth',3);
end;
xlabel('Frame #'); ylabel('PSNR_Y [dB]');

figure,grid,hold on;
for i=1:5,
    plot(sum(akiyo.bits{i}/1000)/length(akiyo.bits{i})*30,meanpsnr(akiyo.psnr{i}.y),fcolor1{i},'MarkerSize',8,'LineWidth',2.2);
end;
legend(legend_str);
xlabel('kbps'); ylabel('PSNR_Y [dB]');

