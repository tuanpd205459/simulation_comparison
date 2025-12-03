
function otf = psf2otf_test(psf,outSize)

if  ~all(psf(:)==0)
   
   psfSize = size(psf);
   % Pad the PSF to outSize
   padSize = outSize - psfSize;
%    psf_PD = zeros(outSize,'gpuArray');
   psf_PD = gpuArray(padarray(psf, padSize, 'post'));

   % Circularly shift otf so that the "center" of the PSF is at the
   % (1,1) element of the array.
   psf_PD    = gpuArray(circshift(psf_PD,-floor(psfSize/2)));

   % Compute the OTF
   otf = gpuArray(fftn(psf_PD));

   % Estimate the rough number of operations involved in the 
   % computation of the FFT.
   nElem = gpuArray(prod(psfSize));
   nOps  = 0;
   for k=1:ndims(psf_PD)
      nffts = nElem/psfSize(k);
      nOps  = nOps + psfSize(k)*log2(psfSize(k))*nffts; 
   end

   % Discard the imaginary part of the psf if it's within roundoff error.
   if max(abs(imag(otf(:))))/max(abs(otf(:))) <= nOps*eps
      otf = gpuArray(real(otf));
   end
else
   otf = zeros(outSize,'gpuArray');
end
end
