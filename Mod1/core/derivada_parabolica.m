% se ajusta una parábola en tres puntos y se toma la derivada en el medio,
% excepto en los extremos, en donde se usan diferencias simples 
% TODO: corregir en los extremos usando parábolas o splines
function D = derivada_parabolica(X, Y)

    assert(length(X) == size(Y,2));
    D = zeros(size(Y));
   
    D(:,1) = diff2(X(1), X(2), X(3), Y(:,1), Y(:,2), Y(:,3), X(1));
    for c=2:length(X)-1
        D(:,c) = diff2(X(c-1), X(c), X(c+1), Y(:,c-1), Y(:,c), Y(:,c+1), X(c));
    end
    D(:,end) = diff2(X(end-2), X(end-1), X(end), Y(:,end-2), Y(:,end-1), Y(:,end), X(end));
    
%     for r=1:size(Y,1)
%         D(r,1) = diff2(X(1), X(2), X(3), Y(r,1), Y(r,2), Y(r,3), X(1));
%         for c=2:length(X)-1
%             D(r,c) = diff2(X(c-1), X(c),X(c+1), Y(r,c-1), Y(r,c), Y(r,c+1), X(c));
%         end
%         D(r,end) = diff2(X(end-2), X(end-1), X(end), Y(r,end-2), Y(r,end-1), Y(r,end), X(end));
%     end
    
end

function d = diff2(x1, x2, x3, y1, y2, y3, v)
    dn  = (x1-x2).*(x1-x3).*(x2-x3);
    d12 = y1-y2;
    d23 = y2-y3;
    d31 = y3-y1;
    b = (x1.^2.*d23 + x2.^2.*d31 + x3.^2.*d12)./dn;
    c =-(x1.*d23 + x2.*d31 + x3.*d12)./dn;
    d = b + 2*c .* v;
end


% function d = diff2(x, y, v)
%     dn= (x(1)-x(2)).*(x(1)-x(3)).*(x(2)-x(3));
%     b = (x(1).^2.*(y(2)-y(3)) + x(2).^2.*(y(3)-y(1)) + x(3).^2.*(y(1)-y(2)))./dn;
%     c = (x(1).*(y(3)-y(2)) + x(2).*(y(1)-y(3)) + x(3).*(y(2)-y(1)))./dn;
%     d = b + 2*c .* v;
% end

