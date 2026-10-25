function [nFibres, nY, nZ] = count_fibres_aligned(cfg, yspace, zspace)
%COUNT_FIBRES_ALIGNED Number of fibres built by mesh_fibres_aligned, without building the mesh.
%   nFibres = count_fibres_aligned(cfg)
%   nFibres = count_fibres_aligned(cfg, yspace, zspace)
%   [nFibres, nY, nZ] = count_fibres_aligned(...)
%
%   mesh_fibres_aligned places one fibre at each (y0, z0) with
%   y0 = cfg.ymin:yspace:cfg.ymax and z0 = cfg.zmin:zspace:cfg.zmax, so
%   nFibres = nY * nZ, with nY = numel(cfg.ymin:yspace:cfg.ymax) and
%   nZ = numel(cfg.zmin:zspace:cfg.zmax) (= floor(range/spacing) + 1 up to
%   round-off). The count does not depend on the rotation angle or on the
%   x-limits. The colon operator is used, as in mesh_fibres_aligned, so the
%   result agrees with it exactly; only 1-D vectors of positions are built.
%
%   Inputs
%     cfg    - scaffold struct (from params_reference): ymin, ymax, zmin,
%              zmax, and yspace, zspace (used when yspace, zspace omitted)
%     yspace - spacings in y: scalar or array
%     zspace - spacings in z: array of the same size as yspace, or a
%              scalar (applied to every yspace)
%
%   Outputs (same size as yspace / zspace)
%     nFibres - number of fibres for each pair (yspace(k), zspace(k))
%     nY, nZ  - number of fibre positions in y and in z
%
%   Example (cfg from params_reference: y in [0.5, 49.5], z in [0.5, 29.5]):
%     count_fibres_aligned(cfg)                              % 28 (= 7*4)
%     count_fibres_aligned(cfg, [7.5 5 10], [9.5 9.5 5])     % [28 40 30]

if nargin < 2
    yspace = cfg.yspace;
    zspace = cfg.zspace;
elseif nargin < 3
    error('count_fibres_aligned:MissingInput', 'Give both yspace and zspace, or neither.');
end

% A scalar spacing is applied to every value of the other one.
if isscalar(yspace)
    yspace = repmat(yspace, size(zspace));
elseif isscalar(zspace)
    zspace = repmat(zspace, size(yspace));
end
if ~isequal(size(yspace), size(zspace))
    error('count_fibres_aligned:SizeMismatch', 'yspace and zspace must have the same size.');
end
if any(yspace(:) <= 0) || any(zspace(:) <= 0)
    error('count_fibres_aligned:InvalidSpacing', 'yspace and zspace must be positive.');
end

% Number of positions per direction, exactly as the colon operator in
% mesh_fibres_aligned generates them.
nY = arrayfun(@(d) numel(cfg.ymin:d:cfg.ymax), yspace);
nZ = arrayfun(@(d) numel(cfg.zmin:d:cfg.zmax), zspace);

nFibres = nY .* nZ;

end
