function hash=pmcHash(text)
md=java.security.MessageDigest.getInstance('SHA-256');
md.update(uint8(unicode2native(text,'UTF-8')));
hash=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
