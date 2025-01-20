% NUFFT_IMPLEMENTATION - non-uniform FFT implementation class
%   A nuFFT implementation is a nuFFT of which is defined for a certain
%   trajectory, onto a certain image size, using a certain oversampling
%   ratio, level of gridding accuracy and resampling method (full or no
%   precomputation).
%
%   nufft_implementation methods:
%   nufft_implementation - class constructor
%   forward - forward nuFFT
%   adjoint - adjoint nuFFT
%   write_impfile - write nuFFT implementation to file
%   delete - clear nuFFT implementation
%   mtimes - overloaded multiplication (*) operator for forward/adjoint nuFFT
%   ctranspose - overloaded transpose (') operator for toggling forward/adjoint nuFFT

% Michal Zarrouk, July 2013.

classdef nufft_implementation < handle
    properties %(SetAccess = private, Hidden = true)
        objectHandle; % Handle to the underlying nuFFT_implementation_t instance - a pointer to the C++ nuFFT_implemetation_t instance
        datatype;
        direction = 'forward';
        Nthreads = 1;
        dontdelete = 0;
    end
    methods
        %% Constructor - Create a new nuFFT_implementation
        function this = nufft_implementation(varargin)
            % NUFFT_IMPLEMENTATION - class constructor
            %   Call sequences:
            %     imp = nufft_implementation(datatype, Nthreads, imagesize, maxerr, alpha, resamplingmethod, trajfilename);
            %     imp = nufft_implementation(datatype, Nthreads, imagesize, maxerr, alpha, resamplingmethod, trajectory, sqrtdcf)
            %     imp = nufft_implementation(datatype, Nthreads, impfilename);
            %     imp = nufft_implementation(datatype, pointer_to_nuFFT_implementation_t);
            % 
            %   Arguments:
            %     datatype - either |'double'| or |'float'|, in consistency
            %       with the class of the input trajectory data.
            %     Nthreads - integer number of threads.
            %     imagesize - array of image size in pixels in each dimension.
            %     maxerr - maximum aliasing error.
            %     alpha - grid oversampling ratio.
            %     resamplingmethod - either |'spm'| (full precomputation) or
            %       |'onthefly'| (no precomputation).
            %     trajfilename - name of trajectory data file, which
            %       contains consecutive sets of k-space coordinates and
            %       dcfs (not sqrt of dcfs!) per each k-space sampling
            %       point, written in binary format, in either double or
            %       float.
            %     trajectory - vector of trajectory coordinates, ordered by
            %       sample.
            %     sqrtdcf - vector of square root of sampling point density
            %       compensation factors.
            %     impfilename - name of a previously-written nuFFT
            %       implementation file.
            %     imp - nuFFT implementation object.
            
            if nargin > 0
                this.datatype = varargin{1};
                if nargin == 2
                    pointer_to_nuFFT_implementation_t = varargin{2};
                    this.objectHandle = pointer_to_nuFFT_implementation_t;
                else
                    this.Nthreads = varargin{2};
                    this.objectHandle = nufft_mex(this, 'init', varargin{2:end});
                end
            end
        end
        
        %% Destructor - Destroy the nufft_implementation class instance
        function delete(this)
            % DELETE - class destructor
            %   imp.delete;
            if ~this.dontdelete
                nufft_mex(this, 'delete', this.objectHandle);    
            end
        end

        %% forward
        function kdata = forward(this, Nthreads, image_data)
            % FORWARD - Forward nuFFT
            %   Call sequence:
            %     kdata = imp.forward(Nthreads, image_data);
            %   
            %   Arguments:
            %     Nthreads - integer number of threads.
            %     image_data - gridded data (may be in a multidimensional
            %       array format).
            %     kdata - vector of transformed data onto non-uniform
            %       sampling locations.
			kdata = nufft_mex(this, 'forward', Nthreads, this.objectHandle, image_data);
        end
        
        %% adjoint
        function image_data = adjoint(this, Nthreads, kdata)
            % ADJOINT - Adjoint nuFFT
            %   Call sequence:
            %     kdata = imp.forward(Nthreads, image_data);
            %   
            %   Arguments:
            %     Nthreads - integer number of threads.
            %     kdata - vector of non-uniformly sampled data.
			%     image_data - gridded data in a multidimensional array
			%       format.
            image_data = nufft_mex(this, 'adjoint', Nthreads, this.objectHandle, kdata);
        end
        
        %% create_impfile
        function write_impfile(this, impfilename)
            % WRITE_IMPFILE - Write nuFFT implementation to an implementation file
            %   Call sequence:
            %     imp.write_impfile(impfilename);
            %
            %   Arguments:
            %     impfilename - name of destination nuFFT implementation
            %       file
            nufft_mex(this, 'write_impfile', this.objectHandle, impfilename);
        end
        
        function res = mtimes(this, data)
            % * - Forward/adjoint nuFFT
            %   Call sequence:
            %     outdata = imp * indata;
            %
            %   Arguments:
            %     imp - nuFFT implementation object
            %     indata - data for transformation
            %     outdata - transformed data
            %       (see adjoint and forward methods for data format
            %       description)
            %    
            %   The default operation for this operator would be a forward
            %   nuFFT, unless imp has been previously transposed using ' .
            
            if strcmp(this.direction,'forward')
                res = this.forward(this.Nthreads,data);
            elseif strcmp(this.direction,'adjoint')
                res = this.adjoint(this.Nthreads,data);
            end
        end
        
        function res = ctranspose(this)
            % ' - Toggle forward/adjoint transformation
            %   Call sequence:
            %     imp2 = imp';
            %     outdata = imp'* indata;
            %
            %   Arguments:
            %     imp - nuFFT implementation object
            %     imp2 - nuFFT implementation where the default operation
            %       for the * operator has been toggled between
            %       forward/adjoint.
            %     indata - data for transformation
            %     outdata - transformed data
            %       (see adjoint and forward methods for data format
            %       description)
            %    
            
            res = copy(this);
            if strcmp(res.direction,'forward')
                res.direction = 'adjoint';
            elseif strcmp(res.direction,'adjoint')
                res.direction = 'forward';
            end
        end
        
        
        function res = copy(this)
            res = nufft_implementation;
            a = properties(this);
            for i=1:length(a)
                res.(a{i}) = this.(a{i});
            end
            res.dontdelete = 1;
        end
        
    end
end


function varargout = nufft_mex(this, varargin)
    if strcmp(this.datatype,'double')
        if nargout == 0
            nufft_mex_double(varargin{:});
            varargout = {};
        else
            varargout{:} = nufft_mex_double(varargin{:});
        end
    elseif strcmp(this.datatype,'float')
        if nargout == 0
            nufft_mex_float(varargin{:});
            varargout = {};
        else
            varargout{:} = nufft_mex_float(varargin{:});
        end
    end

end
