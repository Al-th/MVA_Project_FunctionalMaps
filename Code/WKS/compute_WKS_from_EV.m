% This is the implementation of the Wave Kernel Signature described
% in the paper:
% 
%    The Wave Kernel Signature: A Quantum Mechanical Approach To Shape Analysis 
%    M. Aubry, U. Schlickewei, D. Cremers
%    In IEEE International Conference on Computer Vision (ICCV) - Workshop on 
%    Dynamic Shape Capture and Analysis (4DMOD), 2011
% 
% Please refer to the publication above if you use this software. 
% 
% This work is licensed under a Creative Commons
% Attribution-NonCommercial 3.0 Unported License  
% ( http://creativecommons.org/licenses/by-nc/3.0/ )
% 
% The WKS is patented and violation to the license agreement would
% imply a patent infringement.
%
% THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY EXPRESSED OR IMPLIED WARRANTIES
% OF ANY KIND, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THIS SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THIS SOFTWARE.


function WKS = compute_WKS_from_EV(E,PHI,WKS_options)

% function [WKS,E,PHI] = compute_WKS(vertices,faces) compute
%   the Wave Kernel Signature of triangle mesh given by
%   [vertices,faces]
%   
%   INPUT:
%   vertices is (number of vertices) x 3 matrix
%   faces is a (number of faces) x 3 matrix
%   
%   OUTPUT:
%   WKS is the (number of vertices) x 100 WKS matrix 
%   E is the vector of LB eigenvalues (by default of size 100 x 1)
%   PHI is the (number of vertices x 100) matrix of LB eigenfunctions 
%   L is the cotan Laplace-Beltrami matrix
%
%   The main parameter to adjust depending on your task is wks_variance

if exist('WKS_options','var');
    if ~isfield(WKS_options,'display')
        WKS_options.display=1;
    end
    if ~isfield(WKS_options,'n_eigenvalues')
        WKS_options.n_eigenvalues=100;
    end
    if ~isfield(WKS_options,'N')
        WKS_options.N=100;
    end
    if ~isfield(WKS_options,'variance')
        WKS_options.variance=6;
    end
else
        WKS_options.N=100;
        WKS_options.n_eigenvalues=100;
        WKS_options.display=1;
        WKS_options.variance=6;
end
%% parameters

n_eigenvalues=WKS_options.n_eigenvalues; % number of eigenvalues used for computations
% depending on the application, you can use less than 300
N = WKS_options.N; % number of evaluations of WKS
wks_variance = WKS_options.variance; % variance of the WKS gaussian (wih respect to the 
% difference of the two first eigenvalues). For easy or precision tasks 
% (eg. matching with only isometric deformations) you can take it smaller



%% basic quantities

num_vertices = size(PHI,1);
%num_faces = size(faces,1);

if length(E)~=WKS_options.n_eigenvalues
    if length(E)>WKS_options.n_eigenvalues
        E=E(1:WKS_options.n_eigenvalues);
        PHI=PHI(:,1:WKS_options.n_eigenvalues);
    else
    error('input problem: number of EV');
    end
end




if WKS_options.display

%% compute WKS 

fprintf('Computing WKS...');
end
WKS=zeros(num_vertices,N);

log_E=log(max(abs(E),1e-6))';
e=linspace(log_E(2),(max(log_E))/1.02,N);  
sigma=(e(3)-e(2))*wks_variance;

C = zeros(1,N); %weights used for the normalization of f_E
% 
% PHI2=PHI.^2;
% for i = 1:N
%     WKS(:,i) = PHI2*(exp((-(e(i) - log_E).^2) ./ (2*sigma.^2))');
%     %sum(PHI2.*...
%      %   repmat( exp((-(e(i) - log_E).^2) ./ (2*sigma.^2)),num_vertices,1),2);
%     C(i) = sum(exp((-(e(i)-log_E).^2)/(2*sigma.^2)));
% end

WKS = (PHI.^2)*((exp((-(repmat(e',[1 n_eigenvalues]) - repmat(log_E,[N 1])).^2) ./ (2*sigma.^2)))');
C = sum((exp((-(repmat(e',[1 n_eigenvalues]) - repmat(log_E,[N 1])).^2) ./ (2*sigma.^2)))');

% normalize WKS
WKS(:,:) = WKS(:,:)./repmat(C,num_vertices,1);


if WKS_options.display
fprintf('done. \n');
end


