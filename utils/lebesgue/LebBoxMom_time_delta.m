function moments = LebBoxMom_time_delta(d, box, t, YalmipBasis, normalized)
    if(~exist('YalmipBasis','var') || isempty(YalmipBasis))
        YalmipBasis = 0;
    end

    if nargin < 4
        normalized = 0;
    end
    n = size(box,2);
 
    if(YalmipBasis == 1)
%         disp('Generating moments in Yalmip basis')
        dv = monpowers(n+1,d);
    else
%         disp('Generating moments in Gloptipoly basis')
        dv = genPowGlopti(n+1,d);
    end
    dv_x = dv(:, 2:end);
    moments = zeros(size(dv_x,1),1);
    for i = 1:numel(moments)
        moments(i) = prod((box(2,:).^(dv_x(i,:)+1) - box(1,:).^(dv_x(i,:)+1)) ./ (dv_x(i,:)+1));
    end

    dv_t = dv(:, 1);
    moments_t = t.^dv_t;

    moments = moments .* moments_t;

    if normalized
        moments = moments/moments(1);
    end
end