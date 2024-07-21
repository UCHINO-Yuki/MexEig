function [res1,res2,res3]=mexeig(func,A,in2,in3,in4)
% func : 'compile'... Only compiling mex
%        'syev', 'syevd', 'syevr', 'syevx', 'geev', 'geevx',
%        'sygv', 'sygvd', 'sygvx', 'ggev', 'ggev3', 'ggevx'
% A    : single or double matrix
% in2  : single or double positive definite matrix,
%        'matrix', 'vector',
%        'N', 'P', 'S', 'B' (for BALANC for 'geevx' and 'ggevx')
%        or abstol
% in3  : 'matrix', 'vector',
%        'N', 'P', 'S', 'B' (for BALANC for 'geevx' and 'ggevx')
%        or abstol
% in4  : 'matrix', 'vector', or abstol
% The order of the arguments after in2 (including in2) is arbitrary.
% The arguments after in2 (including in2) are optional.
%
% written ... 2024-07-21 ... Yuki UCHINO

res1=0;
res2=0;
res3=0;

if strcmp(func,'compile')
    disp('Compiling ... ')
    p = which('mexeig');
    p = p(1:end-8);
    try
        mex('-R2018a', [p 'private/dsyev.c'], '-lmwlapack', '-silent', '-outdir', [p 'private']);
        mex('-R2018a', [p 'private/dsyevd.c'], '-lmwlapack', '-silent', '-outdir', [p 'private']);
        mex('-R2018a', [p 'private/dsyevx.c'], '-lmwlapack', '-silent', '-outdir', [p 'private']);
        mex('-R2018a', [p 'private/dsyevr.c'], '-lmwlapack', '-silent', '-outdir', [p 'private']);
        mex('-R2018a', [p 'private/dgeev.c'], '-lmwlapack', '-silent', '-outdir', [p 'private']);
        mex('-R2018a', [p 'private/dgeevx.c'], '-lmwlapack', '-silent', '-outdir', [p 'private']);
        mex('-R2018a', [p 'private/dsygv.c'], '-lmwlapack', '-silent', '-outdir', [p 'private']);
        mex('-R2018a', [p 'private/dsygvd.c'], '-lmwlapack', '-silent', '-outdir', [p 'private']);
        mex('-R2018a', [p 'private/dsygvx.c'], '-lmwlapack', '-silent', '-outdir', [p 'private']);
        mex('-R2018a', [p 'private/dggev.c'], '-lmwlapack', '-silent', '-outdir', [p 'private']);
        mex('-R2018a', [p 'private/dggev3.c'], '-lmwlapack', '-silent', '-outdir', [p 'private']);
        mex('-R2018a', [p 'private/dggevx.c'], '-lmwlapack', '-silent', '-outdir', [p 'private']);

        mex('-R2018a', [p 'private/ssyev.c'], '-lmwlapack', '-silent', '-outdir', [p 'private']);
        mex('-R2018a', [p 'private/ssyevd.c'], '-lmwlapack', '-silent', '-outdir', [p 'private']);
        mex('-R2018a', [p 'private/ssyevx.c'], '-lmwlapack', '-silent', '-outdir', [p 'private']);
        mex('-R2018a', [p 'private/ssyevr.c'], '-lmwlapack', '-silent', '-outdir', [p 'private']);
        mex('-R2018a', [p 'private/sgeev.c'], '-lmwlapack', '-silent', '-outdir', [p 'private']);
        mex('-R2018a', [p 'private/sgeevx.c'], '-lmwlapack', '-silent', '-outdir', [p 'private']);
        mex('-R2018a', [p 'private/ssygv.c'], '-lmwlapack', '-silent', '-outdir', [p 'private']);
        mex('-R2018a', [p 'private/ssygvd.c'], '-lmwlapack', '-silent', '-outdir', [p 'private']);
        mex('-R2018a', [p 'private/ssygvx.c'], '-lmwlapack', '-silent', '-outdir', [p 'private']);
        mex('-R2018a', [p 'private/sggev.c'], '-lmwlapack', '-silent', '-outdir', [p 'private']);
        mex('-R2018a', [p 'private/sggev3.c'], '-lmwlapack', '-silent', '-outdir', [p 'private']);
        mex('-R2018a', [p 'private/sggevx.c'], '-lmwlapack', '-silent', '-outdir', [p 'private']);
    catch
        error('error');
    end
    disp('Compilation Successful.')
    return
end

if isa(A,'double')
    switch func
        case 'syev'
            if nargout <= 1
                if nargin == 2
                    res1=dsyev(A);
                else
                    res1=dsyev(A,in2);
                end
            else
                if nargin == 2
                    [res1,res2]=dsyev(A);
                else
                    [res1,res2]=dsyev(A,in2);
                end
            end
        case 'syevd'
            if nargout <= 1
                if nargin == 2
                    res1=dsyevd(A);
                else
                    res1=dsyevd(A,in2);
                end
            else
                if nargin == 2
                    [res1,res2]=dsyevd(A);
                else
                    [res1,res2]=dsyevd(A,in2);
                end
            end
        case 'syevx'
            if nargout <= 1
                if nargin == 2
                    res1=dsyevx(A);
                elseif nargin == 3
                    res1=dsyevx(A,in2);
                else
                    [in2,in3]=swapping(in2,in3);
                    res1=dsyevx(A,in2,in3);
                end
            else
                if nargin == 2
                    [res1,res2]=dsyevx(A);
                elseif nargin == 3
                    [res1,res2]=dsyevx(A,in2);
                else
                    [in2,in3]=swapping(in2,in3);
                    [res1,res2]=dsyevx(A,in2,in3);
                end
            end
        case 'syevr'
            if nargout <= 1
                if nargin == 2
                    res1=dsyevr(A);
                elseif nargin == 3
                    res1=dsyevr(A,in2);
                else
                    [in2,in3]=swapping(in2,in3);
                    res1=dsyevr(A,in2,in3);
                end
            else
                if nargin == 2
                    [res1,res2]=dsyevr(A);
                elseif nargin == 3
                    [res1,res2]=dsyevr(A,in2);
                else
                    [in2,in3]=swapping(in2,in3);
                    [res1,res2]=dsyevr(A,in2,in3);
                end
            end
        case 'geev'
            if nargout <= 1
                if nargin == 2
                    res1=dgeev(A);
                else
                    res1=dgeev(A,in2);
                end
            elseif nargout == 2
                if nargin == 2
                    [res1,res2]=dgeev(A);
                else
                    [res1,res2]=dgeev(A,in2);
                end
            else
                if nargin == 2
                    [res1,res2,res3]=dgeev(A);
                else
                    [res1,res2,res3]=dgeev(A,in2);
                end
            end
        case 'geevx'
            if nargout <= 1
                if nargin == 2
                    res1=dgeev(A);
                elseif nargin == 3
                    res1=dgeev(A,in2);
                else
                    [in2,in3]=swapchar(in2,in3);
                    res1=dgeev(A,in2,in3);
                end
            elseif nargout == 2
                if nargin == 2
                    [res1,res2]=dgeev(A);
                elseif nargin == 3
                    [res1,res2]=dgeev(A,in2);
                else
                    [in2,in3]=swapchar(in2,in3);
                    [res1,res2]=dgeev(A,in2,in3);
                end
            else
                if nargin == 2
                    [res1,res2,res3]=dgeev(A);
                elseif nargin == 3
                    [res1,res2,res3]=dgeev(A,in2);
                else
                    [in2,in3]=swapchar(in2,in3);
                    [res1,res2,res3]=dgeev(A,in2,in3);
                end
            end
        case 'sygv'
            if nargout <= 1
                if nargin == 3
                    res1=dsygv(A,in2);
                else
                    res1=dsygv(A,in2,in3);
                end
            else
                if nargin == 3
                    [res1,res2]=dsygv(A,in2);
                else
                    [res1,res2]=dsygv(A,in2,in3);
                end
            end
        case 'sygvd'
            if nargout <= 1
                if nargin == 3
                    res1=dsygvd(A,in2);
                else
                    res1=dsygvd(A,in2,in3);
                end
            else
                if nargin == 3
                    [res1,res2]=dsygvd(A,in2);
                else
                    [res1,res2]=dsygvd(A,in2,in3);
                end
            end
        case 'sygvx'
            if nargout <= 1
                if nargin == 3
                    res1=dsygvx(A,in2);
                elseif nargin == 4
                    res1=dsygvx(A,in2,in3);
                else
                    [in3,in4]=swapping(in3,in4);
                    res1=dsygvx(A,in2,in3,in4);
                end
            else
                if nargin == 3
                    [res1,res2]=dsygvx(A,in2);
                elseif nargin == 4
                    [res1,res2]=dsygvx(A,in2,in3);
                else
                    [in3,in4]=swapping(in3,in4);
                    [res1,res2]=dsygvx(A,in2,in3,in4);
                end
            end
        case 'ggev'
            if nargout <= 1
                if nargin == 3
                    res1=dggev(A,in2);
                else
                    res1=dggev(A,in2,in3);
                end
            elseif nargout == 2
                if nargin == 3
                    [res1,res2]=dggev(A,in2);
                else
                    [res1,res2]=dggev(A,in2,in3);
                end
            else
                if nargin == 3
                    [res1,res2,res3]=dggev(A,in2);
                else
                    [res1,res2,res3]=dggev(A,in2,in3);
                end
            end
        case 'ggev3'
            if nargout <= 1
                if nargin == 3
                    res1=dggev3(A,in2);
                else
                    res1=dggev3(A,in2,in3);
                end
            elseif nargout == 2
                if nargin == 3
                    [res1,res2]=dggev3(A,in2);
                else
                    [res1,res2]=dggev3(A,in2,in3);
                end
            else
                if nargin == 3
                    [res1,res2,res3]=dggev3(A,in2);
                else
                    [res1,res2,res3]=dggev3(A,in2,in3);
                end
            end
        case 'ggevx'
            if nargout <= 1
                if nargin == 3
                    res1=dggevx(A,in2);
                elseif nargin == 4
                    res1=dggevx(A,in2,in3);
                else
                    [in3,in4]=swapchar(in3,in4);
                    res1=dggevx(A,in2,in3,in4);
                end
            elseif nargout == 2
                if nargin == 3
                    [res1,res2]=dggevx(A,in2);
                elseif nargin == 4
                    [res1,res2]=dggevx(A,in2,in3);
                else
                    [in3,in4]=swapchar(in3,in4);
                    [res1,res2]=dggevx(A,in2,in3,in4);
                end
            else
                if nargin == 3
                    [res1,res2,res3]=dggevx(A,in2);
                elseif nargin == 4
                    [res1,res2,res3]=dggevx(A,in2,in3);
                else
                    [in3,in4]=swapchar(in3,in4);
                    [res1,res2,res3]=dggevx(A,in2,in3,in4);
                end
            end
        otherwise
            error('invalid input');
    end
else
    switch func
        case 'syev'
            if nargout <= 1
                if nargin == 2
                    res1=ssyev(A);
                else
                    res1=ssyev(A,in2);
                end
            else
                if nargin == 2
                    [res1,res2]=ssyev(A);
                else
                    [res1,res2]=ssyev(A,in2);
                end
            end
        case 'syevd'
            if nargout <= 1
                if nargin == 2
                    res1=ssyevd(A);
                else
                    res1=ssyevd(A,in2);
                end
            else
                if nargin == 2
                    [res1,res2]=ssyevd(A);
                else
                    [res1,res2]=ssyevd(A,in2);
                end
            end
        case 'syevx'
            if nargout <= 1
                if nargin == 2
                    res1=ssyevx(A);
                elseif nargin == 3
                    res1=ssyevx(A,in2);
                else
                    [in2,in3]=swapping(in2,in3);
                    res1=ssyevx(A,in2,in3);
                end
            else
                if nargin == 2
                    [res1,res2]=ssyevx(A);
                elseif nargin == 3
                    [res1,res2]=ssyevx(A,in2);
                else
                    [in2,in3]=swapping(in2,in3);
                    [res1,res2]=ssyevx(A,in2,in3);
                end
            end
        case 'syevr'
            if nargout <= 1
                if nargin == 2
                    res1=ssyevr(A);
                elseif nargin == 3
                    res1=ssyevr(A,in2);
                else
                    [in2,in3]=swapping(in2,in3);
                    res1=ssyevr(A,in2,in3);
                end
            else
                if nargin == 2
                    [res1,res2]=ssyevr(A);
                elseif nargin == 3
                    [res1,res2]=ssyevr(A,in2);
                else
                    [in2,in3]=swapping(in2,in3);
                    [res1,res2]=ssyevr(A,in2,in3);
                end
            end
        case 'geev'
            if nargout <= 1
                if nargin == 2
                    res1=sgeev(A);
                else
                    res1=sgeev(A,in2);
                end
            elseif nargout == 2
                if nargin == 2
                    [res1,res2]=sgeev(A);
                else
                    [res1,res2]=sgeev(A,in2);
                end
            else
                if nargin == 2
                    [res1,res2,res3]=sgeev(A);
                else
                    [res1,res2,res3]=sgeev(A,in2);
                end
            end
        case 'geevx'
            if nargout <= 1
                if nargin == 2
                    res1=dgeev(A);
                elseif nargin == 3
                    res1=dgeev(A,in2);
                else
                    [in2,in3]=swapchar(in2,in3);
                    res1=dgeev(A,in2,in3);
                end
            elseif nargout == 2
                if nargin == 2
                    [res1,res2]=dgeev(A);
                elseif nargin == 3
                    [res1,res2]=dgeev(A,in2);
                else
                    [in2,in3]=swapchar(in2,in3);
                    [res1,res2]=dgeev(A,in2,in3);
                end
            else
                if nargin == 2
                    [res1,res2,res3]=dgeev(A);
                elseif nargin == 3
                    [res1,res2,res3]=dgeev(A,in2);
                else
                    [in2,in3]=swapchar(in2,in3);
                    [res1,res2,res3]=dgeev(A,in2,in3);
                end
            end
        case 'sygv'
            if nargout <= 1
                if nargin == 3
                    res1=ssygv(A,in2);
                else
                    res1=ssygv(A,in2,in3);
                end
            else
                if nargin == 3
                    [res1,res2]=ssygv(A,in2);
                else
                    [res1,res2]=ssygv(A,in2,in3);
                end
            end
        case 'sygvd'
            if nargout <= 1
                if nargin == 3
                    res1=ssygvd(A,in2);
                else
                    res1=ssygvd(A,in2,in3);
                end
            else
                if nargin == 3
                    [res1,res2]=ssygvd(A,in2);
                else
                    [res1,res2]=ssygvd(A,in2,in3);
                end
            end
        case 'sygvx'
            if nargout <= 1
                if nargin == 3
                    res1=ssygvx(A,in2);
                elseif nargin == 4
                    res1=ssygvx(A,in2,in3);
                else
                    [in3,in4]=swapping(in3,in4);
                    res1=ssygvx(A,in2,in3,in4);
                end
            else
                if nargin == 3
                    [res1,res2]=ssygvx(A,in2);
                elseif nargin == 4
                    [res1,res2]=ssygvx(A,in2,in3);
                else
                    [in3,in4]=swapping(in3,in4);
                    [res1,res2]=ssygvx(A,in2,in3,in4);
                end
            end
        case 'ggev'
            if nargout <= 1
                if nargin == 3
                    res1=sggev(A,in2);
                else
                    res1=sggev(A,in2,in3);
                end
            elseif nargout == 2
                if nargin == 3
                    [res1,res2]=sggev(A,in2);
                else
                    [res1,res2]=sggev(A,in2,in3);
                end
            else
                if nargin == 3
                    [res1,res2,res3]=sggev(A,in2);
                else
                    [res1,res2,res3]=sggev(A,in2,in3);
                end
            end
        case 'ggev3'
            if nargout <= 1
                if nargin == 3
                    res1=sggev3(A,in2);
                else
                    res1=sggev3(A,in2,in3);
                end
            elseif nargout == 2
                if nargin == 3
                    [res1,res2]=sggev3(A,in2);
                else
                    [res1,res2]=sggev3(A,in2,in3);
                end
            else
                if nargin == 3
                    [res1,res2,res3]=sggev3(A,in2);
                else
                    [res1,res2,res3]=sggev3(A,in2,in3);
                end
            end
        case 'ggevx'
            if nargout <= 1
                if nargin == 3
                    res1=sggevx(A,in2);
                elseif nargin == 4
                    res1=sggevx(A,in2,in3);
                else
                    [in3,in4]=swapchar(in3,in4);
                    res1=sggevx(A,in2,in3,in4);
                end
            elseif nargout == 2
                if nargin == 3
                    [res1,res2]=sggevx(A,in2);
                elseif nargin == 4
                    [res1,res2]=sggevx(A,in2,in3);
                else
                    [in3,in4]=swapchar(in3,in4);
                    [res1,res2]=sggevx(A,in2,in3,in4);
                end
            else
                if nargin == 3
                    [res1,res2,res3]=sggevx(A,in2);
                elseif nargin == 4
                    [res1,res2,res3]=sggevx(A,in2,in3);
                else
                    [in3,in4]=swapchar(in3,in4);
                    [res1,res2,res3]=sggevx(A,in2,in3,in4);
                end
            end
        otherwise
            error('invalid input');
    end
end
end

%% swapping
function [c,d]=swapping(a,b)
if isa(a,'char')
    c=b;
    d=a;
else
    c=a;
    d=b;
end
end

%% swapchar2
function [c,d]=swapchar(a,b)
if length(a)>2
    c=b;
    d=a;
else
    c=a;
    d=b;
end
end
