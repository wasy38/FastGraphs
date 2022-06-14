using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization.Formatters.Binary;
using System.Text;
using System.Threading.Tasks;

namespace FastGraphs
{
    static class ExtensionMethods
    {
        public static T DeepCopy<T>(this T self) 
        {
            if (!typeof(T).IsSerializable)
                throw new ArgumentException("Type must be serializable");

            if (ReferenceEquals(self, null))
                return default(T);

            var formatter = new BinaryFormatter();
            using (var stream = new MemoryStream())
            {
                formatter.Serialize(stream, self);
                stream.Seek(0, SeekOrigin.Begin);
                return (T)formatter.Deserialize(stream);
            }
        }
    }
}
