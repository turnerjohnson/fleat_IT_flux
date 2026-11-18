function [im, sz] = capture_frame(fh, dpi, sz)
%CAPTURE_FRAME  One movie frame, rendered offscreen at a locked size.
%
%   IM = CAPTURE_FRAME(FH, DPI) renders figure FH to an RGB array at DPI dots
%   per inch and returns it with both dimensions forced even.
%
%   [IM, SZ] = CAPTURE_FRAME(FH, DPI, SZ) additionally checks the frame
%   against SZ, the [rows cols] returned for the first frame, and errors if
%   they differ. Pass SZ = [] on the first frame to establish it.
%
%   This replaces GETFRAME, which reads the figure back OFF THE SCREEN. That
%   readback is at the mercy of things that have nothing to do with the plot:
%   whether another window overlaps the figure, whether the display has gone
%   to sleep, Retina backing-scale rounding, and which monitor the window
%   happens to sit on. Any of those can change the returned size by a pixel
%   or two partway through a render.
%
%   A pixel is enough to ruin the movie. VideoWriter's H.264 encoder fixes
%   its plane strides on the FIRST frame it is given. Hand it a later frame
%   one pixel wider and every row is read at the wrong offset, so each row
%   shifts a little further than the one above it -- the picture shears into
%   diagonal bands. The chroma planes are subsampled, so they shear by a
%   different amount again and the colour separates from the luminance and
%   smears. The encoder does not complain; it just writes garbage. Hence the
%   explicit size check here: a mismatch should stop the render, not survive
%   into the file.
%
%   PRINT -RGBIMAGE renders into an offscreen buffer instead, so none of the
%   screen-state problems apply and the size follows only from the figure's
%   PaperPosition and DPI. The even-dimension crop is belt and braces: H.264
%   subsamples chroma 2x2 and wants even dimensions anyway.
%
%   Note that the caller still needs a plain DRAWNOW before this. DRAWNOW
%   LIMITRATE is free to skip the update entirely, which would capture the
%   previous frame's state.
%
%   See also GETFRAME, PRINT, VIDEOWRITER.

if nargin < 2 || isempty(dpi), dpi = 150; end
if nargin < 3, sz = []; end

im = print(fh, '-RGBImage', sprintf('-r%d', round(dpi)));

% H.264 wants even dimensions; trim rather than pad so nothing is invented.
im = im(1:2*floor(size(im,1)/2), 1:2*floor(size(im,2)/2), :);

now_sz = [size(im,1) size(im,2)];
if isempty(sz)
    sz = now_sz;
elseif ~isequal(now_sz, sz)
    error('capture_frame:sizeChanged', ...
          ['Frame size changed mid-render: [%d %d], expected [%d %d]. ' ...
           'Writing this would shear the movie. Check that nothing is ' ...
           'resizing the figure between frames.'], ...
          now_sz(1), now_sz(2), sz(1), sz(2));
end
